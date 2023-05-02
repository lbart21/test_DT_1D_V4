from Algorithms.DT_1D_V4.post_processing.generate_plots \
    import  GenerateTransientInterfacePropertyPlots, \
            GenerateSinglePlots, \
            GenerateWaterfallPlots, \
            GenerateThrustPlot, \
            GenerateTransientCellPropertyPlots,\
            SingleTemporalPlotsWithMultipleYAxes,\
            SingleSpatialPlotsWithMultipleYAxes

#GenerateTransientCellPropertyPlots(cell_file_name = "Sim23DataAtCellID0.txt", \
                                   #plot_vars = ["p", "T", "rho", "u", "massf", "massf_N2", "massf_N"])

#GenerateSinglePlots(data_file = "Sim5DataAt0.050000033ForComponentRocket.txt", \
                    #plot_vars = ["p", "T", "vel_x", "rho", "massf"])

GenerateWaterfallPlots(data_files = ["Sim3DataAt0.000000000ForComponentrocket.txt", \
                                     "Sim3DataAt0.001000000ForComponentrocket.txt", \
                                    "Sim3DataAt0.003000000ForComponentrocket.txt", \
                                     "Sim3DataAt0.005000000ForComponentrocket.txt"], \
                        plot_vars = ["p", "T", "vel_x", "rho", "Ma", "massf"])

#SingleSpatialPlotsWithMultipleYAxes(data_file = "Sim3DataAt0.100000200ForComponentRocket.txt", \
                                #plot_vars = [["p"], ["Ma"], ["D"]], visible_axes = [1, 1, 0], plot_number = 1)

#GenerateWaterfallPlots(data_files = ["dataAt0.000000000forComponentPipe.txt", \
                                    #"dataAt0.001000018forComponentPipe.txt", \
                                    #"dataAt0.002000347forComponentPipe.txt", \
                                    #"dataAt0.003000202forComponentPipe.txt"], \
                    #plot_vars = ["p", "rho", "vel_x", "Ma", "T"])

GenerateTransientCellPropertyPlots(cell_file_name = "Sim5DataAtCellID0.txt", \
                                    plot_vars = ["p", "rho", "T", "vel_x", "massf"])

#SingleTemporalPlotsWithMultipleYAxes(data_file = "Sim5DataAtCellID0.txt", plot_vars = [["massf_H2", "massf_O2", "massf_H2O"]], visible_axes = [1, 1], plot_number = 1)

#SingleTemporalPlotsWithMultipleYAxes(data_file = "Sim1DataAtCellID0.txt", plot_vars = [["molef_H2", "molef_O2", "molef_H2O"], ["molef_N2"]], visible_axes = [1, 1], plot_number = 2)

#SingleTemporalPlotsWithMultipleYAxes(data_file = "Sim1DataAtCellID0.txt", plot_vars = [["conc_H2", "conc_O2", "conc_H2O"], ["conc_N2"]], visible_axes = [1, 1], plot_number = 3)

#SingleSpatialPlotsWithMultipleYAxes(data_file = "Sim3DataAt0.005000000ForComponentrocket.txt", \
                                    #plot_vars = [["Ma"], ["D"]], visible_axes = [1, 0], plot_number = 3)

#SingleSpatialPlotsWithMultipleYAxes(data_file = "Sim3DataAt0.005000000ForComponentrocket.txt", \
                                    #plot_vars = [["D"]], visible_axes = [1], plot_number = 4)

#GenerateTransientInterfacePropertyPlots(interface_file_name = "Sim3DataAtInterfaceID0.txt", plot_vars = ["mass_flux"])

#GenerateTransientInterfacePropertyPlots(interface_file_name = "Sim3DataAtInterfaceID1.txt", plot_vars = ["mass_flux"])

#GenerateTransientInterfacePropertyPlots(interface_file_name = "Sim3DataAtInterfaceID2.txt", plot_vars = ["mass_flux"])

#GenerateTransientInterfacePropertyPlots(interface_file_name = "Sim3DataAtInterfaceID150.txt", plot_vars = ["Ma", "mass_flux"])