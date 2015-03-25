#!/usr/bin/env python

from pprint import pformat
from os import path, mkdir
from random import gauss
from execo import configuration
from execo_engine import Engine, ParamSweeper, sweep, slugify, logger
from execo.log import style
from numpy import linspace, array, zeros, absolute
from scipy.integrate import odeint
import matplotlib.pyplot as plt 
from solvers import evol_ratio
from plots import plot_isotopic_state

class BoxModel(Engine):
    """
    Base execo engine that perform computation of flux-box equilibrium
    """


    


class IsotopicBoxModel(BoxModel):
    """
    Execo engine that compuet the evolution of an isotopic ratio in
    a flux-box system
    """
    def __init__(self, params=None, boxes=None, flux=None, partcoeff=None):
        """ """
        super(IsotopicBoxModel, self).__init__()
        self.params = params
        self.Boxes = boxes
        self.Flux = flux
        self.Partcoeff = partcoeff

    def run(self):
        """ """
        self.initial_state()
        self._show_initial_state()

    def compute_evolution(self, Delta, func=None, outdir=None):
        """ """
        logger.info(style.log_header('Computing evolution'))
        if func is None:
            func = self.evol_ratio
        if outdir is None:
            outdir = self.result_dir
        Ratio = [(delta / 1e3 + 1e0) * self.standard for delta in Delta]
        Ratio = odeint(func, Ratio, self.time)
        Delta = ((Ratio / self.standard) - 1.0) * 1000
        self.plot_evolution(Delta, outdir=outdir)
        return Delta

    def initial_state(self):
        """ """
        self._show_state('Intial state')
        self._Mass = array([box['mass']
                            for box in self.Boxes.itervalues()])
        self._Flux = array([box.values()
                            for box in self.Flux.values()])
        self._Partcoeff = array([box.values()
                                 for box in self.Partcoeff.values()])
        plot_isotopic_state(self.Boxes.keys(), [b['delta']
                            for b in self.Boxes.values()],
                            self.Flux)

    def _show_state(self, message=None):
        """ """
        if message:
            logger.info(style.user2(message))

#        logger.info('%s \n%s\n%s %s\n%s %s',
#                    style.log_header(message),
#                    ''.join([style.emph(box.rjust(10))
#                             for box in self.Boxes.keys()]),
#                    style.object_repr('Delta'.ljust(10)),
#                    ' '.join(['%f'.ljust(8) % round(box['delta'], 3)
#                             for box in self.Boxes.values()]),
#                    style.object_repr('Mass'.ljust(10)),
#                    ' '.join(['%e' % box['mass']
#                             for box in self.Boxes.values()]))
        logger.info('Boxes\n%s', pformat(self.Boxes))
#        logger.info('Partition coefficients\n%s', pformat(self.Partcoeff))
#        logger.info('Flux\n%s', pformat(self.Flux))
        
        
        
    
#class BoxModel(Engine):
#    """
#    This is the main engine that is used to created custom isotopic models
#    """
#    def __init__(self):
#        """ Add options for the number of measures, migration bandwidth, number of nodes
#        walltime, env_file or env_name, stress, and clusters and initialize the engine """
#        super(BoxModel, self).__init__()
#        self.init_plots()
#        logger.setLevel('INFO')
#        logger.info(style('\n\n                 Welcome to the human isotopic Box Model\n', 'log_header'))
#            
#    def initial_state(self, outdir = None):
#        """ Convert the dict given from parameters to Numpy array """
#        logger.info(style('Initial boxes configuration\n', 'log_header')+
#                    ''.ljust(8)+''.join( [ style(box.rjust(10), 'emph') for box in self.Boxes.iterkeys() ])+
#                    style('\n'+'Delta'.ljust(8), 'object_repr')+''.join( [ str(box['Delta']).rjust(10) for box in self.Boxes.itervalues()])+
#                    style('\n'+'Mass'.ljust(8), 'object_repr')+''.join( [ str(box['Mass']).rjust(10) for box in self.Boxes.itervalues() ])
#                    ) 
#        if outdir is None:
#            outdir = self.result_dir+'/'
#        self.plot_state(self.Boxes.keys(), array( [ box['Delta'] for box in self.Boxes.itervalues()]), 
#                        name = '_initial', outdir = outdir)
#
#        self._Mass = array( [ box['Mass'] for box in self.Boxes.itervalues() ] )
#        self._Flux = array( [ box.values() for box in self.Flux.values() ])
#        self._Partcoeff = array( [ box.values() for box in self.Partcoeff.values() ])
#        
#        f = open(outdir+'Delta.initial', 'w')
#        for box, value in self.Boxes.iteritems():
#            f.write(box+' '+str(value['Delta'])+'\n')
#        f.close()
#        return [ box['Delta'] for box in self.Boxes.itervalues() ]
#        
#    
#    def compute_evolution(self, Delta, func = None, outdir = None):
#        """ """
#        logger.info(style('Computing evolution', 'log_header'))
#        if func is None:
#            func = self.evol_ratio
#        if outdir is None:
#            outdir = self.result_dir
#        Ratio = [ ( delta/1e3+1e0)*self.standard for delta in Delta ] 
#        Ratio = odeint(func, Ratio, self.time)
#        Delta = ((Ratio/self.standard)-1.0)*1000;        
#        self.plot_evolution(Delta, outdir = outdir)
#        return Delta
#    
#    def evol_ratio(self, ratio, t):
#        """ The evolution function that can be used for isotopic ratio evolution"""
#        rationew = zeros(ratio.size)
#        for ii in range(ratio.size):
#            outflux = 0;
#            influx = 0;
#            for jj in range(ratio.size):
#                outflux = outflux + self._Flux[ii][jj]/self._Mass[ii]*self._Partcoeff[ii][jj]*ratio[ii]
#                influx = influx + self._Flux[jj][ii]/self._Mass[ii]*self._Partcoeff[jj][ii]*ratio[jj]
#            rationew[ii]= influx - outflux
#        return rationew;
#    
#    def final_state(self, Delta_final, outdir = None):
#        """ """
#        if outdir is None:
#            outdir = self.result_dir
#        logger.info(style('Final boxes state\n', 'log_header')+
#                    ''.ljust(8)+''.join( [ style(box.rjust(10), 'emph') for box in self.Boxes.iterkeys() ])+
#                    style('\n'+'Delta'.ljust(8), 'object_repr')+''.join( [ str(round(delta, 7)).rjust(10) for delta in Delta_final if absolute(delta) < 1000]))
#        self.plot_state(self.Boxes.keys(), Delta_final, name = '_final', outdir = outdir)
#        f = open(outdir+'Delta.final', 'w')
#        for box in self.Boxes.iterkeys():
#            idx = self.Boxes.keys().index(box)
#            f.write(box+' '+str(Delta_final[idx])+'\n')
#        f.close()
#        
#        
#    def init_plots(self):
#        """ Define the colors and shape of the model boxes"""
#        self.plots_conf = {
#            "diet":   {'color': "#BBFFB5", 'shape': "rectangle"},
#            "plasma": {'color': "#D6DE42", 'shape': "ellipse"},
#            "RBC":    {'color': "#FF3A25", 'shape': "ellipse"},
#            "liver":  {'color': "#93BB8F", 'shape': "ellipse"}, 
#            "urine":  {'color': "#FBFF93", 'shape': "rectangle"},
#            "feces":  {'color': "#BF9285", 'shape': "rectangle"},
#            "menses": {'color': "#FF0600", 'shape': "rectangle"},
#            "kidney": {'color': "#620A00", 'shape': "ellipse"},
#            "muscle": {'color': "#FF2933", 'shape': "ellipse"},
#            "skin":   {'color': "#FFE486", 'shape': "ellipse"} ,
#            "bone":   {'color': "#C9C985", 'shape': "ellipse"}  
#            }
#        self.color_chars = '0123456789ABCDEF'
#        
#    
#    def plot_state(self, boxes, deltas, name = '', outdir = None):
#        """ Make a graph of a given state """
#        graph = Dot(graph_type='digraph', fontname="Verdana", size="10, 5", fixedsize= True)
#        i_box = 0
#        for box in boxes:            
#            textcolor = 'white' if sum( [ self.color_chars.index(col) for col in self.plots_conf[box]['color'].split('#')[1] ] ) < 35 else 'black' 
#            node_box = Node(box, style="filled", label = '<<font POINT-SIZE="10" color="'+textcolor+'">'+box+'<br/> '+
#                            "%.7f" % round(deltas[i_box], 7)+'</font>>',
#                fillcolor = self.plots_conf[box]['color'], shape = self.plots_conf[box]['shape'])
#            i_box += 1
#            graph.add_node(node_box)
#            
#        for box_from, boxes_to in self.Flux.iteritems():
#            for box_to, flux in boxes_to.iteritems():
#                if flux !=0:
#                    if flux > 0:
#                        edge = Edge(box_from, box_to,  label = '<<font POINT-SIZE="10">'+str(flux)+'</font>>')
#                    elif flux < 0:
#                        edge = Edge(box_to, box_from,  label = '<<font POINT-SIZE="10">'+str(flux)+'</font>>')                
#                    graph.add_edge(edge)
#        
#        if outdir is None:
#            outdir = self.result_dir
#            
#        outfile = outdir+'/state'+name+'.png'
#        graph.write_png(outfile)
#        logger.info('State has been saved to '+style(outfile, 'emph'))
#        
#    def plot_evolution(self, Delta, outdir = None):
#        """ Draw a graph of the boxes evolution through years"""
#        plt.figure()
#        i_box = 0
#        for box in self.Boxes:
#            # remove deriving boxes 
#            if absolute(Delta[0,i_box]-Delta[-1,i_box]) < 1000:
#                plt.plot(self.time/365., Delta[:,i_box], label=box,
#                         color = self.plots_conf[box]['color'])
#            i_box += 1
#        plt.legend()
#        plt.xlabel(r"Years")
#        plt.ylabel(self.delta_name)
#        if outdir is None:
#            outdir = self.result_dir
#        outfile = outdir+'/evolution.png'
#        plt.savefig(outfile)
#        logger.info('Evolution has been saved to '+style(outfile, 'emph'))    
#        
#configuration['color_styles']['comb'] = 'on_cyan', 'bold'
#
#
#class IsotopicBoxModel(Engine):
#    """This is the main engine that is used to created custom isotopic models
#    """
#    def __init__(self):
#        """Initialize the execo engine"""
#        super(IsotopicBoxModel, self).__init__()
#        self.init_plots()
#        logger.info(style.log_header('\n\n                 Welcome to the ' +
#                                     'human isotopic Box Model\n'))
#        logger.debug(pformat(self.__dict__))
#
#    def initial_state(self, outdir=None):
#        """ Convert the dict given from parameters to Numpy array """
#        logger.info(style.log_header('Initial boxes configuration\n') +
#                    ''.ljust(8) +
#                    ''.join([style.emph(box.rjust(10))
#                             for box in self.Boxes.iterkeys()]) +
#                    style.object_repr('\n' + 'Delta'.ljust(8)) +
#                    ''.join([str(box['Delta']).rjust(10)
#                             for box in self.Boxes.itervalues()]) +
#                    style.object_repr('\n' + 'Mass'.ljust(8)) +
#                    ''.join([str(box['Mass']).rjust(10)
#                             for box in self.Boxes.itervalues()])
#                    )
#        if outdir is None:
#            outdir = self.result_dir + '/'
#        self.plot_state(self.Boxes.keys(),
#                        array([box['Delta']
#                               for box in self.Boxes.itervalues()]),
#                        name='_initial', outdir=outdir)
#
#        self._Mass = array([box['Mass']
#                            for box in self.Boxes.itervalues()])
#        self._Flux = array([box.values()
#                            for box in self.Flux.values()])
#        self._Partcoeff = array([box.values()
#                                 for box in self.Partcoeff.values()])
#
#        f = open(outdir + '/Delta.initial', 'w')
#        for box, value in self.Boxes.iteritems():
#            f.write(box + ' ' + str(value['Delta']) + '\n')
#        f.close()
#        return [box['Delta'] for box in self.Boxes.itervalues()]
#
#    def compute_evolution(self, Delta, func=None, outdir=None):
#        """ """
#        logger.info(style.log_header('Computing evolution'))
#        if func is None:
#            func = self.evol_ratio
#        if outdir is None:
#            outdir = self.result_dir
#        Ratio = [(delta / 1e3 + 1e0) * self.standard for delta in Delta]
#        Ratio = odeint(func, Ratio, self.time)
#        Delta = ((Ratio / self.standard) - 1.0) * 1000
#        self.plot_evolution(Delta, outdir=outdir)
#        return Delta
#
#    def evol_ratio(self, ratio, t):
#        """ The evolution function that is used for isotopic ratio evolution"""
#        rationew = zeros(ratio.size)
#        for ii in range(ratio.size):
#            outflux = 0
#            influx = 0
#            for jj in range(ratio.size):
#                outflux = outflux + self._Flux[ii][jj] / self._Mass[ii] * \
#                    self._Partcoeff[ii][jj] * ratio[ii]
#                influx = influx + self._Flux[jj][ii] / self._Mass[ii] * \
#                    self._Partcoeff[jj][ii] * ratio[jj]
#            rationew[ii] = influx - outflux
#        return rationew
#
#    def final_state(self, Delta_final, outdir=None):
#        """ """
#        if outdir is None:
#            outdir = self.result_dir
#        logger.info(style.log_header('Final boxes state\n',) + ''.ljust(8) +
#                    ''.join([style.emph(box.rjust(10))
#                            for box in self.Boxes.iterkeys()]) +
#                    style.objec_repr('\n' + 'Delta'.ljust(8)) +
#                    ''.join([str(round(delta, 7)).rjust(10)
#                            for delta in Delta_final
#                            if absolute(delta) < 1000]))
#        self.plot_state(self.Boxes.keys(), Delta_final,
#                        name='_final', outdir=outdir)
#        f = open(outdir + '/Delta.final', 'w')
#        for box in self.Boxes.iterkeys():
#            idx = self.Boxes.keys().index(box)
#            f.write(box + ' ' + str(Delta_final[idx]) + '\n')
#        f.close()
#
#    def init_plots(self):
#        """ Define the colors and shape of the model boxes"""
#        self.plots_conf = {
#            "diet":   {'color': "#BBFFB5", 'shape': "s"},
#            "plasma": {'color': "#D6DE42", 'shape': "o"},
#            "RBC":    {'color': "#FF3A25", 'shape': "o"},
#            "liver":  {'color': "#93BB8F", 'shape': "o"},
#            "urine":  {'color': "#FBFF93", 'shape': "s"},
#            "feces":  {'color': "#BF9285", 'shape': "s"},
#            "menses": {'color': "#FF0600", 'shape': "s"},
#            "kidney": {'color': "#620A00", 'shape': "o"},
#            "muscle": {'color': "#FF2933", 'shape': "o"},
#            "skin":   {'color': "#FFE486", 'shape': "o"},
#            "bone":   {'color': "#C9C985", 'shape': "o"}
#            }
#        self.color_chars = '0123456789ABCDEF'
#
#    def plot_state(self, boxes, deltas, name='', outdir=None):
#        """ Make a graph of a given state """
#        gr = nx.MultiDiGraph()
#        i_box = 0
#        for box in boxes:
#            textcolor = 'white' if sum([self.color_chars.index(col)
#                                        for col in self.plots_conf[box]['color'].split('#')[1]]) < 35 \
#                                else 'black'
#            if box not in gr.nodes():
#                gr.add_node(box, attrib={'delta': round(deltas[i_box], 7)})
#            i_box += 1
#
#        for box_from, boxes_to in self.Flux.iteritems():
#            for box_to, flux in boxes_to.iteritems():
#                if flux != 0:
#                    if flux > 0:
#                        #print box_from, box_to
#                        gr.add_edge(box_from, box_to,
#                                    attrib={'flux': flux})
#                    elif flux < 0:
#                        #print box_to, box_from
#                        gr.add_edge(box_to, box_from,
#                                    attrib={'flux': flux})
#        pos = nx.graphviz_layout(gr, prog='neato')
#
#        for p in gr.nodes():
#            nx.draw_networkx_nodes(gr, pos, nodelist=[p],
#                                   node_color=self.plots_conf[p]['color'],
#                                   node_shape=self.plots_conf[p]['shape'])
##        for f, t, att in gr.edges(data=True):
##            print f, t, att['flux']
##        pos = nx.spring_layout(gr)
##        nx.draw(gr, pos)
##        nx.draw_networkx_labels(gr, pos)
#        plt.axis('off')
#        if outdir is None:
#            outdir = self.result_dir
#
#        outfile = outdir + '/state' + name + '.png'
#
#        plt.savefig(outfile)
#        plt.close()
#        logger.info('State has been saved to ' + style.emph(outfile))
#
#    def plot_evolution(self, Delta, outdir=None):
#        """ Draw a graph of the boxes evolution through years"""
#        fig = plt.figure()
#        i_box = 0
#        for box in self.Boxes:
#            # remove deriving boxes
#            if absolute(Delta[0, i_box] - Delta[-1, i_box]) < 1000:
#                plt.plot(self.time / 365., Delta[:, i_box], label=box,
#                         color=self.plots_conf[box]['color'])
#            i_box += 1
#        plt.legend()
#        plt.xlabel(r"Years")
#        plt.ylabel(self.delta_name)
#        if outdir is None:
#            outdir = self.result_dir
#        outfile = outdir + '/evolution.png'
#        plt.savefig(outfile)
#        logger.info('Evolution has been saved to ' + style.emph(outfile))
#        fig.clf()
#        plt.close()
#        gc.collect()
