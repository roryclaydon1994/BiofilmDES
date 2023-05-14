"""
    Plot all cells on ax using patch collections for speed
"""
# Standard
from matplotlib.collections import PatchCollection, LineCollection

# Custom
from RodShapedBacteria import RodShapedBacterium
from ChainingRodShapedBacteria import ChainingRodShapedBacterium

# Third party
from descartes import PolygonPatch
from shapely.geometry import Point


def add_polygon_to_plot(ax, cell, fc, ec):
    cell.addElementToPlot(ax, colour=fc, ec=ec)


def addAllCellsToPlot(cells, ax, ax_rng, alpha=0.8, show_id=False, ec='k'):
    """
        Try to optimise adding all cells to an axis
    """

    def cellPlotDict(cell):
        plot_radius = 0.5 * cell.diameter * RodShapedBacterium.sus_vis_radius_factor
        cdict = {
            'polygon': cell.getPolygon(plot_radius),
            'fc': cell.colour,
            'ec': ec,
            'alpha': alpha,
            'lw': 2 / ax_rng
        }
        return cdict

    try:
        patches = [cellPatch(**cellPlotDict(cell)) for cell in cells]
        ax.add_collection(PatchCollection(patches, match_original=True, zorder=1))
    except Exception as e:
        for cell in cells:
            cell.colour = RodShapedBacterium.colour_fn(cell.theta)
            add_polygon_to_plot(ax=ax, cell=cell, fc=cell.colour, ec='w')


def addSpringLinks(cells, ax, ax_rng, show_id=False):
    def findSprings(cell):
        for id, params in cell.springs.items():
            cell2 = cell.cell_hash[id]
            p1 = cell.rcm + params[0] * cell.ori * cell.length * 0.5
            p2 = cell2.rcm + params[1] * cell2.ori * cell2.length * 0.5
            ov = np.array(params[2:])
            # This needs to be implemented once correct surface point selected
            # cv=p2-p1
            p1 += cell.sus_vis_radius_factor * cell.radius * ov
            p2 -= cell2.sus_vis_radius_factor * cell2.radius * ov
            return np.array([p1[:2], p2[:2]])

    lines = [findSprings(cell) for cell in cells if cell.springs]
    anchors = [Point(pnt).buffer(0.1) for line in lines for pnt in line]
    patches = [PolygonPatch(a, fc=None, ec='k', alpha=0.8, lw=5 / ax_rng)
               for a in anchors
               ]
    connected_fraction = len([True for cell in cells if cell.springs]) / len(cells)
    print(f"{connected_fraction=}")
    ax.add_collection(
        LineCollection(lines, lw=5 / ax_rng, colors='k')
    )
    ax.add_collection(
        PatchCollection(patches, match_original=True)
    )
    return connected_fraction


def addChainLinks(cells, ax, ax_rng, show_id=False, colour='k'):
    def findChains(cell):
        cell2 = cell.cell_hash[int(cell.upper_link)]
        links = ChainingRodShapedBacterium.getLinkVector(cell, cell2, False)
        return [lnk[:, :-1] for lnk in links]

    lines = [findChains(cell) for cell in cells if cell.upper_link != "None"]
    lines = [ll for line in lines for ll in line]
    anchors = [Point(pnt).buffer(0.1) for line in lines for pnt in line]
    patches = [PolygonPatch(a, fc=colour, ec=colour, alpha=0.8, lw=5 / ax_rng)
               for a in anchors
               ]
    # connected_fraction = len([True for cell in cells if cell.upper_link])/len(cells)
    # print(f"{connected_fraction=}")
    ax.add_collection(
        LineCollection(lines, lw=20 / ax_rng, colors=colour)
    )
    ax.add_collection(
        PatchCollection(patches, match_original=True)
    )
