"""Ripped from pyfish on 09 27 2022 to include custom color mapping
https://bitbucket.org/schwarzlab/pyfish/src/main/"""
"""Library to create Fish (Muller) plots."""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import numpy as np
import pandas as pd
from matplotlib import cm
from scipy.ndimage import gaussian_filter


def _stackplot(x, *args, ax=None, colors=None, labels=(), **kwargs):
    """Draw a stacked area plot.

    Taken largely from the matplotlib implementation except that the keywords `edgecolor` and
    `where` are provided to `fill_between`
    """
    y = np.row_stack(args)
    labels = iter(labels)

    if ax is None:
        ax = plt.gca()

    if colors is not None:
        ax.set_prop_cycle(color=colors)

    # Assume data passed has not been 'stacked', so stack it here.
    # We'll need a float buffer for the upcoming calculations.
    stack = np.cumsum(y, axis=0, dtype=np.promote_types(y.dtype, np.float32))

    # Color between x = 0 and the first array.
    color = ax._get_lines.get_next_color()
    coll = ax.fill_between(x, 0, stack[0, :], where=(stack[0, :] != 0), interpolate=True,
                           facecolor=color, edgecolor=color, label=next(labels, None),
                           **kwargs)
    coll.sticky_edges.y[:] = [0]
    r = [coll]

    # Color between array i-1 and array i
    for i in range(len(y) - 1):
        color = ax._get_lines.get_next_color()
        r.append(ax.fill_between(x, stack[i, :], stack[i + 1, :],
                                 where=(stack[i, :] != stack[i + 1, :]), facecolor=color,
                                 edgecolor=color, label=next(labels, None), interpolate=True,
                                 **kwargs))

    return r


def _create_ordering(cur_tree, cur_clone):
    """Create index for population DataFrame.

    Recursively traverses the parent tree.
    Build a list such that children are listed between two instances of parent.
    """
    res = []
    if cur_clone in cur_tree.keys():
        res += [cur_clone]
        for child in cur_tree[cur_clone]:
            res += _create_ordering(cur_tree, child)
        res += [cur_clone]
    else:
        return [cur_clone]
    return res


def _build_tree(parent_df):
    """Create a dict-based tree where for each parent there is a list of children."""
    parents = parent_df["ParentId"].unique()
    children = parent_df["ChildId"].unique()
    root_list = np.setdiff1d(parents, children)
    if len(root_list) != 1:
        raise ValueError("Failed to determine root. "
                        "There must be exactly one node with no parent in the parent tree.")
    root_id = root_list[0]
    ids = np.concatenate([[root_id], children])

    tree = {cid: list(children) for cid in ids if len(
        children := parent_df.loc[parent_df["ParentId"] == cid, "ChildId"]) > 0}
    if -1 not in tree:
        tree[-1] = [root_id]

    return tree, ids, root_id


def process_data(pops_df, parent_df,
                 first_step=None, last_step=None, interpolation=0, absolute=False, smooth=None, seed=0,
                 cmap_name="rainbow", color_by=None, colorDict=None):
    """Load data required for plotting.

    Args:
        pops_df (DataFrame): Populations across steps (Id: +int, Step: +int, Pop: +int)
        parent_df (DataFrame): Parent-child relationships (ParentID: +int, ChildID: +int)
        first_step (int, optional): First step to plot. Defaults to None.
        last_step (int, optional): Last step to plot. Defaults to None.
        interpolation (int, optional): Order of interpolation. Defaults to 0.
        absolute (bool, optional): Does not normalize data if set True. Defaults to False.
        smooth (int, optional): Window for Gaussian smoothing. Defaults to None.
        seed (int, optional): Seed used for coloring. Defaults to 0.
        cmap_name (str, optional): Matplotlib colormap. Defaults to None.

    Returns:
        pd.DataFrame: DataFrame with population information
        pd.Series: Series used for x-axis
        list: List of colors
        int: Maximum population. None if not absolute
    """
    min_step = pops_df["Step"].min()
    first_step = max(first_step, min_step) if first_step else min_step

    max_step = pops_df["Step"].max()
    last_step = min(last_step, max_step) if last_step else max_step

    pops_df = pops_df.loc[(pops_df["Step"] >= first_step) & (pops_df["Step"] <= last_step)]

    # populations dataframe processing
    pops_table = pops_df.set_index(['Id', 'Step'])['Pop'].unstack('Step')

    if interpolation > 0:
        if pops_table.shape[1] > 50:
            print("WARNING: interpolation is not recommened for large data")
        pops_table[first_step].fillna(0, inplace=True)
        pops_table[last_step].fillna(0, inplace=True)

        if (~pops_table.isna()).sum(axis=1).min() - 1 < interpolation:
            raise ValueError(f"For interpolation order {interpolation}, the iput data has not "
                             f"enough datapoints (at least {interpolation + 1} per sample)")

        larger_pops_table = pd.DataFrame(index=pops_table.index,
                                         columns=np.arange(first_step, last_step - 0.1, 0.1))
        larger_pops_table[pops_table.columns.astype(float)] = pops_table
        larger_pops_table = larger_pops_table.astype(float)

        linear_interpolation = larger_pops_table.interpolate(axis=1, method='linear')
        pops_table_interpolate = larger_pops_table.interpolate(axis=1, method='spline',
                                                               order=interpolation)
        pops_table_interpolate[linear_interpolation <= 0] = 0

        pops_table = pops_table_interpolate

    else:
        pops_table = pops_table.fillna(0)
        if not absolute and (pops_table == 0).all(axis=0).any():
            raise ValueError("If --interpolation is not set and --absolute is not set you cannot "
                             "have missing Steps")

    steps = pops_table.columns

    if smooth:
        pops_table = pd.DataFrame(gaussian_filter(pops_table, (0, smooth)),
                                  index=pops_table.index, columns=pops_table.columns)

    if absolute:
        pops_sums = pops_table.sum(axis=0)
        pop_max = pops_sums.max()
        pops_rest = pop_max - pops_sums
        pops_table.loc[-1] = pops_rest
    else:
        pop_max = None

    # Build parental relationship
    tree, ids, root_id = _build_tree(parent_df)
    samples = pops_df['Id'].unique()
    ordering = _create_ordering(tree, -1 if absolute else root_id)
    ordering = [x for x in ordering if x in samples or x == -1 and absolute]

    pops_stack = pops_table.loc[ordering]
    val, count = np.unique(pops_stack.index, return_counts=True)
    doubles = val[count > 1]
    pops_stack.loc[doubles] = pops_stack.loc[doubles] / 2

    pops_sum = pops_stack.sum(axis=0)
    pops_sum[pops_sum == 0] = 1
    pops_stack = pops_stack / pops_sum

    colors = _create_colors(ids, root_id, ordering, seed, cmap_name, pops_df, color_by=color_by, colorDict=colorDict)
    return pops_stack, steps, colors, pop_max


def setup_figure(width=1920, height=1080, absolute=False, axisLabelFontsize=28, tickLabelFontsize=24, xlabel='Time point'):
    """Create figure with specified height and width."""
    dpi = (width + height) // 20
    plt.figure(figsize=(width // dpi, height // dpi), dpi=dpi)
    plt.xlabel(xlabel, fontsize=axisLabelFontsize)
    plt.ylabel('Population size' if absolute else 'Proportion', fontsize=axisLabelFontsize)
    plt.tick_params(axis='both', labelsize=tickLabelFontsize)
    return


def fish_plot(pops_stack, steps, colors=None, pop_max=None, ax=None, xticksMapping=None):
    """Plot the actual fish plot."""
    _stackplot(steps, pops_stack, colors=colors, ax=ax)

    if ax is None:
        ax = plt.gca()

    first_step = pops_stack.columns[0]
    last_step = pops_stack.columns[-1]
    ax.set_xlim(first_step, last_step)
    ax.set_ylim(0, 1)
    label_text = np.round(np.abs(np.arange(-0.5, 0.6, 0.1)), 1)
    ax.set_yticks(np.arange(0, 1.1, step=0.1))
    ax.set_yticklabels((label_text * pop_max).astype(int) if pop_max is not None else label_text)
    ax.set_xticks(np.arange(first_step, last_step + 1, step=(last_step - first_step) / 10).astype(int))
    plt.draw()
    labels = ax.get_xticklabels()
    if not xticksMapping is None:
        labels = [xticksMapping[str(tick.get_text()).replace('−', '-')] for tick in labels]
    ax.set_xticklabels(labels)
    return

    
def _create_colors(ids, root_id, ordering, seed, cmap_name, pops_df, color_by=None, colorDict=None):
    
    try:
        cmap = cm.get_cmap(cmap_name)
    except ValueError:
        print("WARNING: colormap not recognized, setting to default")
        cmap = cm.get_cmap("rainbow")

    if color_by is not None:
        # color by separate column in pops_df
        assert color_by in pops_df.columns, "'color_by' has to be a column in pops_df"

        min_value = pops_df[color_by].min()
        unique_colors = cmap(np.linspace(0, 1, pops_df[color_by].max() - min_value + 1))
        unique_colors = {value: unique_colors[value - min_value] for value in np.sort(pops_df[color_by].unique())}
        
        colors = pd.DataFrame(np.zeros((len(ids), 4)), index=ids)
        for cur_id in ids:
            colors.loc[cur_id] = unique_colors[pops_df.loc[pops_df['Id']==cur_id, color_by].iloc[0]]

    else:
        ## color by id
        #np.random.seed(seed)
        #colors = np.array(cmap(np.linspace(0, 1, len(ids))))
        #np.random.shuffle(colors)
        #colors = pd.DataFrame(colors, index=ids)
        #colors.loc[-1] = np.ones(4)
        #colors.loc[root_id] = .5 * np.ones(4)
        
        if colorDict is None:
            raise ValueError
        
        colors = pd.DataFrame([matplotlib.colors.to_rgba(colorDict[str(id)]) for id in ids], index=ids)
        colors.loc[root_id] = .5 * np.ones(4)

    colors = colors.loc[ordering].values

    return colors