#!/usr/bin/env python3
#############################################
# This script is to generate color code for the consistency of microbial taxonomy
#############################################

from load import *
class color:
    """
    A class for color generation for visualization
    Random colors
    """
    def __init__(self,i=10):
        self.i = i
    def default_c(self):
        #https://simpledits.com/color-palette-guide-with-pantone-colors-for-spring-summer-2020-nyfw-with-hex-cmyk-and-rgb-values/
        #blue, yellow, green, chive, red
        # faded denim, orange peel, mosaic blue, sunlight, coral pink
        #cinamon stick, grape compote, lark, navy blazer, brilliant white
        #ash
        cc = cycler('color',
                    ['#34558b', '#ffaf12', '#4ec5a5', '#565d47','#d13b40',
                     '#798fa8', '#fd823e', '#117893', '#f0daa4', '#eaac9d',
                     '#a2553a', '#72617d','#b49c73','#3b3d4b','#eff0f1',
                     '#a09d9c',])
        # cc = cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22',
        #              '#17becf'])
        c_list = []
        for j, c in zip(range(self.i), cc()):
            c_list.append(c['color'])
        return c_list

# get the color code
def get_color(all_dict,c_list, name='genus'):
        """
        Assign color to genus/species level

        input
        all_dict[phase_taxonomy] = nxd feature

        output
        row_colors: dict for clustermap color code
        tax_to_color: dict for Raman tree color code
        """
        taxonomy = set()
        for strain in all_dict.keys():
            taxonomy.add(
                strain.split('_')[1].split(' ')[0]
            )
        tax_to_color = dict(zip(taxonomy, c_list))
        row_colors = [
            (strain, tax_to_color[strain.split('_')[1].split(' ')[0]]) \
            for strain in all_dict.keys() if 'Exp' not in strain
        ]
        row_colors = pd.DataFrame(row_colors, columns=['index', name])
        row_colors = row_colors.set_index(row_colors.iloc[:, 0]).drop('index', axis=1)
        return row_colors, tax_to_color

def plot_color(tax_to_color,tax=None):
    # generate color code figure
    temp_color = pd.DataFrame.from_dict(tax_to_color,orient='index')

    temp_color = pd.concat([temp_color,pd.DataFrame(range(temp_color.shape[0]),\
                                                    index=temp_color.index)],axis=1,ignore_index=True)
    temp_color.columns = ['color','value']
    cmap = colors.ListedColormap(temp_color.iloc[:,0].values)
    boundaries = [temp_color.iloc[:,1].values]
    # norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)
    fig,ax = plt.subplots(figsize=(0.5,temp_color.shape[0]/3))
    sns.heatmap(temp_color.iloc[:,[1]],cmap=cmap,ax=ax,
               edgecolor='white',
               linewidth=5,
                cbar=False
               )#,norm=norm)
    ax.set_yticklabels(ax.get_yticklabels(),fontsize='large',weight='bold',style='italic')
    ax.set_xticklabels('')
    rc('text', usetex=False)
    if tax!=None:
        ax.set_title('%s level'%tax,fontsize='x-large',weight='bold')
    return fig
