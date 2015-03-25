import networkx as nx
import matplotlib.pyplot as plt
from execo import logger
from execo.log import style

plots_conf = {"DIET":   {'color': "#BBFFB5", 'shape': "s"},
              "PLASMA": {'color': "#D6DE42", 'shape': "o"},
              "RBC":    {'color': "#FF3A25", 'shape': "o"},
              "LIVER":  {'color': "#93BB8F", 'shape': "o"},
              "URINE":  {'color': "#FBFF93", 'shape': "s"},
              "FECES":  {'color': "#BF9285", 'shape': "s"},
              "MENSES": {'color': "#FF0600", 'shape': "s"},
              "KIDNEY": {'color': "#620A00", 'shape': "o"},
              "MUSCLE": {'color': "#FF2933", 'shape': "o"},
              "SKIN":   {'color': "#FFE486", 'shape': "o"},
              "BONE":   {'color': "#C9C985", 'shape': "o"}}
color_chars = '0123456789ABCDEF'


def plot_isotopic_state(boxes, deltas, flux, name='', outdir=None):
    """ Make a graph of a given state """
    logger.info('Drawing state')
    gr = nx.DiGraph()
    i_box = 0
    for box in boxes:
        textcolor = 'white' if sum([color_chars.index(col)
                                    for col in plots_conf[box]['color'].split('#')[1]]) < 35 \
                            else 'black'
        if box not in gr.nodes():
            gr.add_node(box, attrib={'delta': round(deltas[i_box], 7)})
        i_box += 1

    for box_from, boxes_to in flux.iteritems():
        for box_to, fx in boxes_to.iteritems():
            if fx != 0:
                if fx > 0:
                    print box_from, box_to
                    gr.add_edge(box_from, box_to,
                                attrib={'flux': flux})
                elif flux < 0:
                    print box_to, box_from
                    gr.add_edge(box_to, box_from,
                                attrib={'flux': flux})
    pos = nx.graphviz_layout(gr, prog='dot')

    for p in gr.nodes():
        nx.draw_networkx_nodes(gr, pos, nodelist=[p],
                               node_color=plots_conf[p]['color'],
                               node_shape=plots_conf[p]['shape'])
        
        nx.draw_networkx_edges(gr, pos, nodelist=[p],
                               node_color=plots_conf[p]['color'],
                               node_shape=plots_conf[p]['shape'])
#        for f, t, att in gr.edges(data=True):
#            print f, t, att['flux']
#        pos = nx.spring_layout(gr)
#        nx.draw(gr, pos)
        nx.draw_networkx_labels(gr, pos, arrows=True)
    plt.axis('off')
    
#    if outdir is None:
#        outdir = self.result_dir

#    outfile = outdir + '/state' + name + '.png'
    outfile = "test.png"

    plt.savefig(outfile)
    plt.close()
    logger.info('State has been saved to ' + style.emph(outfile))

def plot_evolution(self, Delta, outdir=None):
    """ Draw a graph of the boxes evolution through years"""
    fig = plt.figure()
    i_box = 0
    for box in self.Boxes:
        # remove deriving boxes
        if absolute(Delta[0, i_box] - Delta[-1, i_box]) < 1000:
            plt.plot(self.time / 365., Delta[:, i_box], label=box,
                     color=self.plots_conf[box]['color'])
        i_box += 1
    plt.legend()
    plt.xlabel(r"Years")
    plt.ylabel(self.delta_name)
    if outdir is None:
        outdir = self.result_dir
    outfile = outdir + '/evolution.png'
    plt.savefig(outfile)
    logger.info('Evolution has been saved to ' + style.emph(outfile))
    fig.clf()
    plt.close()
        