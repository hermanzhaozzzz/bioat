import math
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
def plot_colortable(colors, *, ncols=4, sort_colors=True, labels=None):

    cell_width = 212
    cell_height = 22
    swatch_width = 48
    margin = 12

    names = list(colors)

    if not labels:
        labels = names

    n = len(names)
    nrows = math.ceil(n / ncols)

    width = cell_width * 4 + 2 * margin
    height = cell_height * nrows + 2 * margin
    dpi = 72

    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin / width, margin / height,
                        (width - margin) / width, (height - margin) / height)
    ax.set_xlim(0, cell_width * 4)
    ax.set_ylim(cell_height * (nrows - 0.5), -cell_height / 2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()

    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height

        swatch_start_x = cell_width * col
        text_pos_x = cell_width * col + swatch_width + 7

        ax.text(text_pos_x, y, labels[i], fontsize=14,
                horizontalalignment='left',
                verticalalignment='center')

        ax.add_patch(
            Rectangle(xy=(swatch_start_x, y - 9), width=swatch_width,
                      height=18, facecolor=name, edgecolor='0.7')
        )

    return fig

if __name__ == '__main__':
    # %%% test func plot_colortable
    colors = ['#64C1E8',
              '#80CED7',
              '#63C7B2',
              '#8E6C88',
              '#CA61C3',
              '#FF958C',
              '#883677']
    plot_colortable(colors, ncols=1, labels=[1, 2, 3, 4, 5, 6, 7])
    plt.show()