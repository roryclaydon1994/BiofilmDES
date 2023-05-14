"""
    Chaining RodShapedBacterium element class

    To do:
        Remove recursion from the functions for finding chain head and tail.

    Author: Rory Claydon
"""

# Standard libraries
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg
import copy
from RodShapedBacteria import RodShapedBacterium
import networkx as nx
from shapely.geometry import LineString

class ChainingRodShapedBacterium(RodShapedBacterium):
    """
        Control visualisation of chaining susceptible in a biofilm

    """

    bending_moduli=1
    # cell_hash = {}

    def __init__(self,cell_id,length,radius,
                      pos_x,pos_y,pos_z,
                      ori_x,ori_y,ori_z,
                      upper_link=None,lower_link=None,
                      virtual_center_x=0.0,
                      virtual_center_y=0.0,
                      virtual_center_z=0.0,
                      force_x=0.0,
                      force_y=0.0,
                      force_z=0.0,
                      torque_x=0.0,
                      torque_y=0.0,
                      torque_z=0.0,
                      neighbours=None):
        super().__init__(
            cell_id           = cell_id,
            length            = length,
            radius            = radius,
            pos_x             = pos_x,
            pos_y             = pos_y,
            pos_z             = pos_z,
            ori_x             = ori_x,
            ori_y             = ori_y,
            ori_z             = ori_z,
            virtual_center_x  = virtual_center_x,
            virtual_center_y  = virtual_center_y,
            virtual_center_z  = virtual_center_z,
            force_x           = force_x,
            force_y           = force_y,
            force_z           = force_z,
            torque_x          = torque_x,
            torque_y          = torque_y,
            torque_z          = torque_z,
            neighbours        = neighbours
        )

        self.spring_energy=0
        self.upper_link = upper_link
        self.lower_link = lower_link

    @staticmethod
    def getLinkVector(cell_1,cell_2,primary_only=True):
        la=cell_1.getUpperLinkAnchors()
        ua=cell_2.getLowerLinkAnchors()
        # p1 = cell_1.rcm + 0.5*( cell_1.length ) * cell_1.ori
        # p2 = cell_2.rcm - 0.5*( cell_2.length ) * cell_2.ori
        t0=la['b']-la['a']
        t1=ua['c']-la['b']
        t2=ua['d']-ua['c']
        if primary_only:
            return np.array([la['b'],ua['c']])
        else:
            return (np.array([la['a'],la['b']]),
                    np.array([la['b'],ua['c']]),
                    np.array([ua['c'],ua['d']]))

    @staticmethod
    def addLink(ax,cell_1,cell_2,direction="upper",colour="k"):
        link=getLinkVector(cell_1,cell_2)
        primary_top=link[0]
        secondary_bottom=link[1]
        link_n = primary_top-secondary_bottom
        link_n /= np.linalg.norm( link_n )
        perp_vector = np.array([link_n[1],-link_n[0]])
        # print(self.ori)
        points = [primary_top[0:-1]     - 0.05*perp_vector,
                  primary_top[0:-1]     + 0.05*perp_vector,
                  secondary_bottom[0:-1]+ 0.05*perp_vector,
                  secondary_bottom[0:-1]- 0.05*perp_vector]
        line = plt.Polygon(points, closed=True, facecolor=colour,
                                   edgecolor=None)
        ax.add_patch(line)

        for x,y,_ in [primary_top,secondary_bottom]:
            circle = plt.Circle((x,y), radius=0.1, fc=colour)
            ax.add_patch(circle)

        # points = [primary_top[0:-1],
        #           secondary_bottom[0:-1]]
        # line = plt.Polygon(points,closed=True,facecolor=colour,
        #                    edgecolor=colour,linewidth=3
        #                   )
        # ax.add_patch(line)

        # ax.text(cell_1.rcm[0],
        #         cell_1.rcm[1],
        #         f"{cell_1.cell_id}"
        #        )

    @staticmethod
    def addLinksToElementsInPlot(ax,cell_list,colour='k'):
        """
            Draw the links between susceptible on the plot
            Parameters:
                ax: axis handle
                    The axis to draw the particle on
                cell_list: list of ChainingRodShapedBacterium
                    Add the links to the elements in this list
                colour: colour
                    Specify link fill colour for this element
            Returns:
                perp_vector: float array, 3 elements
                    The vector perpendicular to the orientaion
        """

        # Convert to hash list for fast lookup
        lst_to_dist = lambda lst: {
            lst[ii].cell_id: lst[ii] for ii in range(0,len(lst))
        }
        # cell_hash = lst_to_dist(cell_list)

        for cell_id,cell in ChainingRodShapedBacterium.cell_hash.items():
            if cell.upper_link != "None":
                # print( cell_hash[int(cell.upper_link)] )
                ChainingRodShapedBacterium.addLink(
                    ax,cell,self.cell_hash[int(cell.upper_link)],colour
                )

    @staticmethod
    def getChainHead(cell_id,cell_hash_list):
        upper_link_id = cell_hash_list[int(cell_id)].upper_link
        if upper_link_id == "None":
            return int(cell_id)
        else:
            return ChainingRodShapedBacterium.getChainHead(upper_link_id,cell_hash_list)

    @staticmethod
    def getChainTail(cell_id,cell_hash_list):
        lower_link_id = cell_hash_list[int(cell_id)].lower_link
        if lower_link_id == "None":
            return int(cell_id)
        else:
            return ChainingRodShapedBacterium.getChainTail(lower_link_id,cell_hash_list)

    @staticmethod
    def getChains(cells):

        pairs=[(int(cc.cell_id),int(cc.upper_link)) for cc in cells
                if cc.upper_link != "None"]

        G=nx.Graph()
        G.add_nodes_from([int(cc.cell_id) for cc in cells])
        G.add_edges_from(pairs)

        largest_cc  = max(nx.connected_components(G), key=len)
        smallest_cc = min(nx.connected_components(G), key=len)
        print(f"biggest chain {len(largest_cc)} smallest chain {len(smallest_cc)}")

        chains=[]
        # cell_hash={ int(cell.cell_id) : cell for cell in cells }
        for s in nx.connected_components(G):
            chains.append([ ChainingRodShapedBacterium.cell_hash[id]
                            for id in s
                            ])

        # for s in nx.isolates(G):
        #     chains.append([cell_hash[s]])
        # chains = []
        # print(f'\rFound {len(chains)}', end='', flush=True)
        #
        # while len(cell_hash)>0:
        #
        #     # Populate chains one by one
        #     current_chain = []
        #
        #     # Get next available cell
        #     cell_id,cell = list(cell_hash.items())[0]
        #
        #     # Find chain ends
        #     head = ChainingRodShapedBacterium.getChainHead(cell_id,cell_hash)
        #     tail = ChainingRodShapedBacterium.getChainTail(cell_id,cell_hash)
        #
        #     current_cell = copy.deepcopy(cell_hash[head])
        #
        #     # Traverse the chain until we run out of cells
        #     while True:
        #
        #         current_chain.append(current_cell)
        #
        #         if current_cell.lower_link == "None":
        #             chains.append(current_chain)
        #             del cell_hash[current_cell.cell_id]
        #             break
        #         else:
        #             next_cell_id = int(current_cell.lower_link)
        #             del cell_hash[current_cell.cell_id]
        #             current_cell = copy.deepcopy(cell_hash[next_cell_id])
        #
        #     print(f'\rFound {len(chains)} chains', end='', flush=True)
        # print("\n")
        return chains

    @staticmethod
    def colourChains(ax,cell_list,colours,projection="xy",ax_rng=20):
        """
            Find chains, colour them in according to the average of the head
            and tail
            Parameters:
                ax: axis handle
                    The axis to draw the particle on
                cell_list: list of ChainingRodShapedBacterium
                    Add the links to the elements in this list
                colours: list of colours
                    Specify element fill colour for this chain
                projection: str
                    determines onto which plane to project the 3D biofilm
                    options: "xy","xz","yz"
            Returns:
                perp_vector: float array, 3 elements
                    The vector perpendicular to the orientaion
        """

        # Convert to hash list for fast lookup
        # lst_to_dist = lambda lst: {
        #     lst[ii].cell_id: lst[ii] for ii in range(0,len(lst))
        # }
        # cell_hash = copy.deepcopy(lst_to_dist(cell_list))
        tmp_cell_hash = copy.deepcopy(self.cell_hash)

        max_z = max([cell.pos_z for cell in cell_list])
        max_cap = 15
        max_z = max(max_z,max_cap)

        min_z = min([cell.pos_z for cell in cell_list])
        min_z = max(min_z,0)
        # print(f"{min_z},{max_z}")
        # evenly_spaced_interval = np.linspace(0, 1, 1000)
        # colours = [colour_scheme(x) for x in evenly_spaced_interval]


        # colour_index = 0
        while len(tmp_cell_hash)>0:
            cell_id,cell = list(tmp_cell_hash.items())[0]
            # print(cell_id,cell)
            head = ChainingRodShapedBacterium.getChainHead(cell_id,tmp_cell_hash)
            tail = ChainingRodShapedBacterium.getChainTail(cell_id,tmp_cell_hash)
            # print("The head and tail of this chain is {:d} to {:d}".format(head,tail))

            current_cell = copy.deepcopy(tmp_cell_hash[head])
            # print("Adding: ",current_cell,
            #       "upper: ",current_cell.upper_link,
            #       "lower: ",current_cell.lower_link)
            # quit()
            # print(head%len(colours))
            while True:
                # current_cell.addElementToPlot(ax,
                #                               colour=colours[head%len(colours)],
                #                               projection=projection)
                # print(current_cell)
                capped_colour=min(max(current_cell.pos_z,min_z),max_z)-min_z
                colour = 1-np.log(1+capped_colour)/np.log(1+max_z-min_z)
                # print(f"colour {colour}")
                current_cell.addElementToPlot(
                    ax,
                    colour=colours(colour),
                    projection=projection,
                    ax_rng=ax_rng
                )

                if current_cell.lower_link == "None":
                    del tmp_cell_hash[current_cell.cell_id]
                    break
                else:
                    next_cell_id = int(current_cell.lower_link)
                    del tmp_cell_hash[current_cell.cell_id]
                    current_cell = copy.deepcopy(tmp_cell_hash[next_cell_id])


            # colour_index+=1

            # if cell.upper_link != "None":
            #     # print( cell_hash[int(cell.upper_link)] )
            #     ChainingRodShapedBacterium.addLink(
            #         ax,cell,cell_hash[int(cell.upper_link)],colour
            #     )

    def getLowerLinkAnchors(self):
        c=self.rcm - 0.5*self.ori*self.length
        d=c+self.ori*self.diameter
        return { 'c' : c, 'd' : d }

    def getUpperLinkAnchors(self):
        b=self.rcm + 0.5*self.ori*self.length
        a=b-self.ori*self.diameter
        return { 'a' : a, 'b' : b }

    def getHeadSpringEnergy(self,bend_only=False):
        if self.upper_link == "None":
            return None,None

        upper_cell=self.cell_hash[int(self.upper_link)]
        lower_cell=self

        upper_links=upper_cell.getLowerLinkAnchors()
        lower_links=lower_cell.getUpperLinkAnchors()
        spring_pos=0.5*(lower_links['b']+upper_links['c'])

        t0=lower_links['b']-lower_links['a']
        t1=upper_links['c']-lower_links['b']
        t2=upper_links['d']-upper_links['c']
        links=[t0,t1,t2]

        def getLinkEnergy(l1,l2):
            l1_n=np.linalg.norm(l1)
            l2_n=np.linalg.norm(l2)
            ds=0.5*( l1_n+l2_n )
            kb=2*(np.cross(l1,l2))/(l1_n*l2_n + np.dot(l1,l2))
            return np.dot(kb,kb)/ds

        energy_contributions=[
            getLinkEnergy(links[ii],links[ii+1]) for ii in range(len(links)-1)
            ]
        spring_energy=0.5*self.bending_moduli*np.sum(energy_contributions)
        extension=np.linalg.norm(t1)-(upper_cell.radius+lower_cell.radius)
        if bend_only==False:
            spring_energy+=0.5*self.kappa*extension**2
        spring_energy*=0.5 # Account for double counting

        upper_cell.spring_energy+=spring_energy
        lower_cell.spring_energy+=spring_energy

        return spring_energy,spring_pos[:2]

    @staticmethod
    def getSpringEnergies(cells,bend_only=False):
        # ChainingRodShapedBacterium.makeHashList(cells)
        # N=len(cell_list)
        # spring_energies=np.zeros(N)
        # spring_location=np.zeros((N,2))
        #
        # for cc,cell in enumerate(cell_list):
        #     if cell.upper_link == "None":
        #         continue
        #     spring_energies[cc],spring_location[cc]=cell.getHeadSpringEnergy()

        for cc,cell in enumerate(cells):
            if cell.upper_link != "None":
                cell.getHeadSpringEnergy(bend_only=bend_only)
        return [cell.spring_energy for cell in cells]

    @staticmethod
    def computeChainDescriptors(cells):
        """
            Calculate all statistical properties of the chains. These include:
                Average length
                Average persistence length
        """
        chain_stats = {}
        chains = ChainingRodShapedBacterium.getChains(cells)
        chain_stats['Nchains'] = len(chains)
        chain_stats['n_c'] = np.mean([len(cc) for cc in chains])
        (chain_stats['average_length'],
         chain_stats['average_length_std']
         ) = ChainingRodShapedBacterium.computeAverageLength(chains)
        (chain_stats['persistence_length']
         ) = ChainingRodShapedBacterium.computePersistenceLength(chains)
        return chains, chain_stats

    @staticmethod
    def computeAverageLength(chains):
        """
            Find the average and standard deviation of the chain lengths
        """
        chain_lens=[]
        for chain in chains:
            tmp_len = 0 # length of current chain
            for cell in chain:
                tmp_len += cell.length+2*cell.radius
            chain_lens.append(tmp_len)
        chain_lens = np.array(chain_lens)
        return chain_lens.mean(), chain_lens.std()

    @staticmethod
    def getCurvature(chain):
        """
        Find theta**2 along the chain
        """
        thetas=[]
        for ii in range(len(chain)-1):
            costheta=chain[ii].ori.dot(chain[ii+1].ori)
            if costheta>1:
                costheta=1
            elif costheta<-1:
                costheta=-1
            thetas.append(np.arccos(costheta))
        return np.array(thetas)**2

    @staticmethod
    def computePersistenceLength(chains):
        """
            https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6419224/
        """
        bond_lens=0
        cos_corrs=[]
        num_cells=0
        for chain in chains:
            tmp_len=0
            tmp_corr=0
            cell_0 = chain[0]
            t_0 = cell_0.ori
            l_0 = cell_0.length+2*cell_0.radius
            num_cells+=len(chain)
            for cell in chain:
                l_i = cell.length+2*cell.radius
                tmp_corr += t_0.dot(cell.ori)*l_0*l_i
                bond_lens += l_i
            cos_corrs.append(tmp_corr)
        avg_bond = bond_lens/num_cells
        print(f"{num_cells=}")
        p_len = np.mean(cos_corrs) / avg_bond
        return p_len

    @staticmethod
    def computeEnergyGrid(cell_list,dr=0.5,bend_only=False):

        ChainingRodShapedBacterium.getSpringEnergies(cell_list,bend_only)

        for cell in cell_list:
            cell.energy=cell.spring_energy

        energy_field, energy_grid = RodShapedBacterium.computeEnergyGrid(
            cell_list,
            dr=dr
            )

        return energy_field, energy_grid

    @staticmethod
    def getAllSpringGons():

        if not ChainingRodShapedBacterium.cell_hash:
            print("Error! Please create cell hash list first")

        springs=[]
        for cell_id,cell in ChainingRodShapedBacterium.cell_hash.items():
            if cell.upper_link != "None":
                # print( cell_hash[int(cell.upper_link)] )
                p1=cell.getUpperLinkAnchors()['b']
                p2=ChainingRodShapedBacterium.cell_hash[int(cell.upper_link)].getLowerLinkAnchors()['c']
                springs.append(LineString([p1[:2],p2[:2]]).buffer(0.1))

        return springs

    def getHeadSpringGon(self,hash):
        """
        Only get springs which are in the local hash
        """
        try:
            p1=self.getUpperLinkAnchors()['b']
            p2=hash[int(self.upper_link)].getLowerLinkAnchors()['c']
            return LineString([p1[:2],p2[:2]]).buffer(0.1)
        except Exception as e:
            return None
