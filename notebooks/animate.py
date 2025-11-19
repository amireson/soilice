import matplotlib.pyplot as pl
from matplotlib.animation import FuncAnimation
from IPython.display import HTML

def animate(
    data, z, t,
    animation_configs,
    figsize=(8, 4),
    interval=100,
    blit=True,
    suptitle="Time = {t:.1f} days"
):
    """
    Create a flexible animation with multiple subplots and lines.

    Parameters
    ----------
    data : dict
        Dictionary of named 2D arrays (shape: [time, depth]).
        Keys correspond to variable names used in subplot configs.
    z : array-like
        Depth or y-axis values.
    t : array-like
        Time array (one entry per frame).
    subplot_configs : list of dict
        Each dict describes a subplot, e.g.:
        {
            "vars": ["thetaI3", "thetaL3"],  # variables to plot
            "colors": ["c", "b"],
            "labels": ["Ice", "Liquid"],
            "xlabel": "Water content",
            "ylabel": "Depth (m)",
            "legend": True,
            "xlim": (0, 0.3),
            "ylim": (2, 0)
        }
    figsize : tuple
        Figure size in inches.
    interval : int
        Animation frame interval (ms).
    blit : bool
        Whether to use blitting for faster rendering.
    suptitle : str
        Title format string; may use {t} placeholder for time.
    """

    n_subplots = len(animation_configs)
    fig, axes = pl.subplots(1, n_subplots, figsize=figsize)
    if n_subplots == 1:
        axes = [axes]

    all_lines = []

    # --- Setup each subplot ---
    for ax, cfg in zip(axes, animation_configs):
        n_vars = len(cfg["vars"])
        colors = cfg.get("colors", [None]*n_vars)
        labels = cfg.get("labels", [None]*n_vars)

        lines = [
            ax.plot([], [], lw=2, color=colors[i], label=labels[i])[0]
            for i in range(n_vars)
        ]
        all_lines.extend(lines)

        ax.set_xlim(*cfg.get("xlim", (0, 1)))
        ax.set_ylim(*cfg.get("ylim", (1, 0)))  # depth downward
        ax.set_xlabel(cfg.get("xlabel", ""))
        ax.set_ylabel(cfg.get("ylabel", ""))

        if cfg.get("legend", False):
            ax.legend(loc=cfg.get("legend_loc", "best"))

        if ax != axes[0]:
            ax.set_yticklabels('')
            
        ax.grid(True)

    # --- Init function ---
    def init():
        for line in all_lines:
            line.set_data([], [])
        return all_lines

    # --- Update function ---
    def update(frame):
        line_index = 0
        for cfg in animation_configs:
            for var in cfg["vars"]:
                all_lines[line_index].set_data(data[var][frame, :], z)
                line_index += 1
        fig.suptitle(suptitle.format(t=t[frame]))
        return all_lines

    ani = FuncAnimation(fig, update, frames=len(t),
                        init_func=init, blit=blit, interval=interval)

    pl.close(fig)
    return ani
