import glob
import matplotlib.pyplot

class PlotInterpolation(object):
    """
    Make plots to support plotting of the PPSolveFactoryTest output
    """
    def __init__(self, file_name_list):
        """
        Initialise the plotter - file_name_list is the list of files that will
        be plotted
        """
        self.file_name_list = file_name_list
        self.data = {}
        self.plots = {}
        self.setup_plots()

    def setup_plots(self):
        """
        Set up empty plots
        """
        self.plots["fit_fig"] = matplotlib.pyplot.figure(0, (12.8, 9.6))
        self.plots["fit_fig"].xlabel = "s"
        self.plots["fit_fig"].ylabel = "value"
        kwd = {"ylim":[-1, 1], "xlim":[-1.0, 20.0]}
        ax = self.plots["fit_fig"].add_subplot(1, 1, 1, **kwd)
        self.plots["residual_fig"] = matplotlib.pyplot.figure(1, (12.8, 9.6))
        self.plots["residual_fig"].xlabel = "s"
        self.plots["residual_fig"].ylabel = "residual"
        kwd = {"yscale":"log", "ylim":[1e-6, 1000], "xlim":[-1.0, 25.0]}
        ax = self.plots["residual_fig"].add_subplot(1, 1, 1, **kwd)

    def get_fit_fig(self):
        matplotlib.pyplot.figure(0)
        return self.plots["fit_fig"]

    def get_residual_fig(self):
        matplotlib.pyplot.figure(1)
        return self.plots["residual_fig"]

    def finish_plots(self):
        self.get_fit_fig()
        matplotlib.pyplot.legend()
        self.get_residual_fig()
        matplotlib.pyplot.legend()
        matplotlib.pyplot.show(block=False)
        self.plots["fit_fig"].savefig("ppsolvefactorytest_fit.png")
        self.plots["residual_fig"].savefig("ppsolvefactorytest_residual.png")

    def load_file(self, file_name):
        fin = open(file_name)
        self.data = {}
        self.data["title"] = fin.readline().rstrip()
        smooth = self.data["title"].split("smooth")[1]
        smooth = smooth.split("_")[1]
        poly = self.data["title"].split("poly")[1]
        poly = poly.split("_")[1]
        self.data["smooth"] = int(smooth)
        self.data["poly"] = int(poly)
        columns = fin.readline()
        columns.rstrip()
        self.data["columns"] = columns.split()
        for col in self.data["columns"]:
            self.data[col] = []
        for line in fin.readlines():
            line.rstrip()
            for i, a_word in enumerate(line.split()):
                col = self.data["columns"][i]
                self.data[col].append(float(a_word))
        self.data["length"] = len(self.data[col])
        self.get_s_data()
        self.get_fit_errors()

    def get_s_data(self):
        s_list = [0.]
        for i in range(self.data["length"]):
            x, y, z = self.data["pos_x"][i], self.data["pos_y"][i], self.data["pos_z"][i]
            if i == 0:
                print("Stepping from", x, y, z, end=' ')
            if i > 0:
                s = ((x-x_old)**2+(y-y_old)**2+(z-z_old)**2)**0.5
                s_list.append(s+s_list[-1])
            x_old, y_old, z_old = x, y, z
        print("to", x, y, z)
        self.data["s_pos"] = s_list

    def get_fit_errors(self):
        for key in "x", "y", "z":
            true = self.data["true_"+key]
            fit = self.data["fit_"+key]
            error = [abs(fit[i]-true[i]) for i in range(self.data["length"])]
            self.data["residual_"+key] = error

    def get_title_suffix(self):
        return "t: "+ str(self.data["smooth"])+ " f: "+str(self.data["poly"])

    def do_fit_plots(self):
        self.get_fit_fig()
        title = "Fitted - "+self.get_title_suffix()
        abs_list = self.data["s_pos"]
        ord_list = self.data["fit_x"]
        matplotlib.pyplot.plot(abs_list, ord_list,
                                    color=self.color(), marker=self.marker(), 
                                    fillstyle='none', label=title)

    def color(self):
        i = self.data['smooth']
        c_list = ['red', 'blue', 'green', 'orange', 'magenta']
        color = c_list[i%len(c_list)]
        return color

    def marker(self):
        i = self.data['poly']
        m_list = ['o', 'v', '^', 's', 'x']
        marker = m_list[i%len(m_list)]
        return marker

    def do_truth_plot(self):
        self.get_fit_fig()
        abs_list = self.data["s_pos"]
        ord_list = self.data["true_x"]
        color = self.color()
        matplotlib.pyplot.plot(abs_list, ord_list, label="True x")

    def do_residual_plot(self):
        self.get_residual_fig()
        title = "Residual - "+self.get_title_suffix()
        abs_list = self.data["s_pos"]
        ord_list = self.data["residual_x"]
        matplotlib.pyplot.semilogy(abs_list, ord_list,
                                    color=self.color(), marker=self.marker(), 
                                    fillstyle='none', label=title)

    def main(self):
        for i, file_name in enumerate(sorted(self.file_name_list)):
            self.load_file(file_name)
            self.do_fit_plots()
            if i == len(self.file_name_list)-1:
                self.do_truth_plot()
            self.do_residual_plot()
        self.finish_plots()


def main():
    plotter = PlotInterpolation(glob.glob("PPSolveFactoryTest_*"))
    plotter.main()
    input("Press <CR> to end")

if __name__ == "__main__":
    main()