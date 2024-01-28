# -*- coding: utf-8 -*-
# pylint: disable=E0012, fixme, invalid-name, no-member, R0913, R1716, W0613, W0622, E0611, W0603, R0902, R0903, C0301
"""  Semester 2, Assignment 4

Assignment Tasks: 8

Aim:
    GUI to visualise electric field lines for a collection of point charges.

    To remove a charge:
        Left-click on it and release.

    To add a charge:
        Left-click and hold, move mouse left or right reading charge in status
        message and release. Adds the charge of displayed magnitude at the initial
        mouse press position.

    Drag a charge:
        Left-click on a charge and drag it

Restrictions:
    Do not change anything outside TODO blocks.
    Do not use import.
    Do not add pylint directives.

Guidance:
    Please read the introduction page to this assignment on Canvas.

Author of the template:
    Wolfgang Theis
    School of Physics and Astronomy
    University of Birmingham
"""

import sys
import numpy as np
from scipy.integrate import odeint
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget
from PyQt5.QtWidgets import QSizePolicy, QVBoxLayout
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas



class Charges:
    ''' A collection of point charges in the xy-plane

    distances, charges, electric fields, and potentials are treated as unitless

    The electric field due to a charge q at distance d is taken to be
    E = q/d**2

    The potential V is given by dV/dr = -E

    '''

    def __init__(self):
        # private attributes
        # self._q is a 1d array holding the value of the charge for the different charges
        # charge n has a charge of self._q[n]
        self._q = np.zeros((0, ))
        # self._pos is a 2d array with xy coordinates of the charges
        # x position of charge n is self._pos[n, 0]
        # y position of charge n is self._pos[n, 1]
        self._pos = np.zeros((0, 2))

    def add_charge(self, q, xy):
        ''' add charge of magnitude q at locations x,y = xy '''
        self._q = np.hstack([self._q, q])
        self._pos = np.vstack([self._pos, np.array([xy])])

    def delete_charge(self, k):
        ''' delete charge k '''
        if 0 <= k < self._q.shape[0]:
            self._q = np.delete(self._q, k)
            self._pos = np.delete(self._pos, k, 0)

    def set_position(self, k, xy):
        ''' set position of charge k to xy '''
        if 0 <= k < self._q.shape[0]:
            self._pos[k, :] = xy

    def get_charges(self):
        ''' provide list of (charge, position) tuples '''
        p = self._pos
        return [(q, p[k, :]) for k, q in enumerate(self._q)]

    def get_closest(self, xy, limit=0.1):
        ''' determine the index of the closest charge

        Parameters
        ----------
        xy: array-like
            x,y pair of position to find the closest charge to.

        limit:
            maximum distance to find a charge in.

        Returns
        -------
        index: int
            index of the closest charge
            If no charge is within the limit then None is returned.
        '''
        # TODO: Assignment Task 1: write function body
        x_pos = self._pos[0:, 0]
        y_pos = self._pos[0:, 1]
        # get all the charges in the array, work out their distances from the clicked position,
        # see if the minimum in the list is within 0.1 of the clicked position
        charge_distances = np.sqrt((xy[0]-x_pos)**2 + (xy[1]-y_pos)**2)
        if np.min(charge_distances) <= limit:
            # if within limit, return the index
            return np.argmin(charge_distances)
        return None
        # End of Task 1; proceed to task 2.

    def scaled_electric_field(self, xy, _):
        ''' calculate suitably scaled electric field vector at position xy.

        The scaling is such that integrating this scaled field along field lines
        Es(lambda) dlambda from lambda=0 to lambda_max transverses roughly a
        modulus of the potential energy difference of lambda_max.

        To avoid numerical instabilities it should be ensured that the return value
        is not divergent by ensuring that the scaling factor remains finite.
        If you had a pure scaling factor of 1/x where x could be zero,
        you would use 1/(x+0.0001) instead]

        Parameters
        ----------
        xy: array-like
            x,y pair of position at which the scaled electric field is requested

        _:
            not used but required place holder for ode integrator

        Returns
        -------
        ef: array-like
            scaled electric field Es = (s*Ex, s*Ey) at position xy

        '''
        # TODO: Assignment Task 2: write function body
        # for all charges, add their electric field vectors, then scale them down
        charge_values = self._q
        x_pos = self._pos[0:, 0]
        y_pos = self._pos[0:, 1]
        distances = np.array((np.sqrt((xy[0]-x_pos)**2 + (xy[1]-y_pos)**2)))
        Ex = ([charge_values]/(distances**3)) * (xy[0]-x_pos)
        Ey = ([charge_values]/(distances**3)) * (xy[1]-y_pos)
        E_total = [np.sum(Ex), np.sum(Ey)]
        E_scaled = E_total/(E_total[0]**2 + E_total[1]**2 + 0.001)
        return E_scaled
        # End of Task 2; proceed to task 3.

    def field_lines(self, nr_of_fieldlines=32, start_radius=0.2, lambda_max=10, points=801):
        ''' calculate the field lines which should include one at pi/4, rather than 0

            Parameters
            ----------

            nr_of_fieldlines: int
                number of field lines originating from each positive charge

            start_radius: float
                radius around each positive charge at which the field lines start

            lambda_max: float
                the maximum value of the parameter for the parametric representation
                of the fieldlines x(lambda), y(lambda), lambda =[0,..., lambda_max]

            points: int
                the number of points in each fieldline

            Returns
            -------

            fieldlines: list of 2-d numpy arrays
                each element of the list is an array of shape (points, 2) with
                fieldlines[k][:, 0] and fieldlines[k][:, 1] providing
                the x and y values, respectively of the k-th fieldline
        '''
        # TODO: Assignment Task 3: write function body
        lambda_array = np.linspace(0, lambda_max, points)
        lines = []
        for q, pos in self.get_charges():
            if q > 0:
                # unpacking the xy coordinates of the centre of the charge
                x, y = pos
                for f in range(0, nr_of_fieldlines):
                    # calculating the angle between each starting point in radians
                    theta = 2 * f * np.pi / nr_of_fieldlines
                    # adding the extra bits that take you round the edge of the
                    # charge in this (nr_of_fieldlines) number of steps
                    x0 = x + np.cos(theta)*start_radius
                    y0 = y + np.sin(theta)*start_radius
                    # integrating along each line and adding it to the list of lines
                    lines.append(odeint(self.scaled_electric_field, [x0, y0], lambda_array))
        return lines
        # End of Task 3; proceed to task 4.


class MyMainWindow(QMainWindow):
    ''' the main window potentially with menus, statusbar, ... '''

    def __init__(self):
        super().__init__()
        self.resize(500, 500)
        self.move(400, 300)
        central_widget = MyCentralWidget(self)
        self.setCentralWidget(central_widget)
        self.setWindowTitle('PyQt widget with matplotlib figure')
        self.statusBar().showMessage('Waiting for mouse move or mouse click')


class MyCentralWidget(QWidget):
    ''' everything in the main area of the main window '''

    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        # define figure canvas
        self.mpl_widget = MyMplWidget(self.main_window)
        # place MplWidget into a vertical box layout
        vbox = QVBoxLayout()
        vbox.addWidget(self.mpl_widget)  # add the figure
        # use the box layout to fill the window
        self.setLayout(vbox)


class MyMplWidget(FigureCanvas):
    ''' both a QWidget and a matplotlib figure '''

    def __init__(self, main_window, parent=None, figsize=(4, 4), dpi=100):
        self.main_window = main_window
        self.fig = plt.figure(figsize=figsize, dpi=dpi)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        self.mouse_pressed_pos = None
        self.closest_k = None
        self.dragging = False
        self.qadd = 0
        self.lines = []
        self.points = []
        self.field_lines_args = None            # used to save the parameters for use in drag_replt
        self.charges = Charges()
        # add some charges to start with an example
        self.charges.add_charge(1, (1, 0))
        self.charges.add_charge(1, (-1, 0))
        self.charges.add_charge(-1, (0, 1))
        self.charges.add_charge(-1, (0, -1))
        self.plot_fieldlines()
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.cid_press = self.fig.canvas.mpl_connect('button_press_event', self.on_mouse_press)
        self.cid_release = self.fig.canvas.mpl_connect('button_release_event', self.on_mouse_release)
        self.cid_move = self.fig.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)

    def plot_fieldlines(self, nr_of_fieldlines=32, start_radius=0.2, lambda_max=10, points=801):
        ''' clear figure and plot fieldlines and charges '''
        self.fig.clf()
        self.ax = self.fig.add_subplot(1, 1, 1)
        self.lines = []         # list of matplotlib lines in the plot showing the fieldlines (for drag_replot)
        self.points = []        # list of matplotlib lines in the plot showing the charges (for drag_replot)
        self.field_lines_args = (nr_of_fieldlines, start_radius, lambda_max, points)
        # TODO: Assignment Task 4: calculate and plot field lines; plot charges; collect lines and points
        # the shape of this is no of charges by nr_of_fieldlines by points by 2
        # so first loop through each charge, each line on the charge, each point(on the line)'s xy coordinates
        field_lines = self.charges.field_lines(nr_of_fieldlines, start_radius, lambda_max, points)
        for line in field_lines:
            line_plot, = self.ax.plot(line[:, 0], line[:, 1], 'k', zorder=0)
            self.lines.append(line_plot)
            for i in self.charges.get_charges():
                if i[0] < 0:
                    # negative charge, so it's red
                    color = 'r'
                else:
                    # positive charge, so it's blue
                    color = 'b'
                self.points.append(self.ax.add_artist(plt.Circle(i[1], 0.2, color=color, zorder=10)))
        # End of Task 4; proceed to task 5.
        self.ax.set_xlabel('$x$')
        self.ax.set_ylabel('$y$')
        self.ax.set_xlim(-5, 5)
        self.ax.set_ylim(-5, 5)
        self.ax.set_aspect('equal')
        self.draw()

    def drag_replot(self):
        ''' replot elements in self.lines and self.point
            after adjusting their xdata and ydata arrays
            to reflect the new position of the dragged charge.
        '''
        self.ax.draw_artist(self.ax.patch)                              # <-- redraw the plotting area of the axis
        # TODO: Assignment Task 5: redraw updated field lines and charges
        # calculate new 8 lines
        field_lines = self.charges.field_lines(8, self.field_lines_args[1], self.field_lines_args[2], self.field_lines_args[3])
        # fieldlines need to be drawn with nr_of_fieldlines = 8 not 32
        # so only plot every 4th line
        for i, line in enumerate(self.lines[::4]):
            line.set_xdata(field_lines[i][:, 0])
            line.set_ydata(field_lines[i][:, 1])
            self.ax.draw_artist(line)
        for point in self.points:
            self.ax.draw_artist(point)
        # End of Task 5; proceed to task 6.
        self.fig.canvas.update()                                        # <-- update the figure
        self.fig.canvas.flush_events()                                  # <-- ensure all draw requests are sent out

    def on_mouse_move(self, event):
        ''' add charge or drag

            self.dragging determines whether a charge is being dragged and is
            updated accordingly.

            If a charge is being added then self.qadd, the value of the charge
            is updated and the current value displayed in the statusbar.

            If a charge is being dragged then its position is updated in the
            self.charges object and the charges and fieldlines are displayed using
            drag_replot().

        '''
        # TODO: Assignment Task 6: write function body
        # if the mouse is on the plot, let the user know where it is
        if event.xdata is not None and event.ydata is not None:
            msg = f'x={event.xdata:.2f}, y={event.ydata:.2f}'
            if self.mouse_pressed_pos is not None: # if the mouse has been clicked,
                # check if a charge was clicked
                if self.closest_k is not None:
                    # if the mouse has clicked on an existing charge,
                    #set dragging to true
                    self.dragging = True
                    # drag replot called
                    self.charges.set_position(self.closest_k, [event.xdata, event.ydata])
                    self.drag_replot()
                    # update the status message with the current mouse xy coords
                    self.main_window.statusBar().showMessage(msg)
                else:
                    # if the mouse hasn't clicked on an existing charge
                    # then a new one is being added, so calculate
                    # the charge of the new one
                    self.qadd = round(event.xdata - self.mouse_pressed_pos[0])
                    # update the status message with the current xy coords and
                    # the charge of the new charge
                    msg_q = f', q={self.qadd:.2f}'
                    self.main_window.statusBar().showMessage(msg + msg_q)
        elif event.xdata is not None and event.ydata is not None:
            # even if the mouse hasnt been clicked, the status message should
            # still show the xy position in the status message
            self.main_window.statusBar().showMessage(msg)
        # End of Task 6; proceed to task 7.

    def on_mouse_press(self, event):
        ''' set self.mouse_pressed_pos and self.closest_k, the
            xy pair for the position the mouse was clicked at
            and the index of of the charge closest to that position, respectively.
        '''
        # TODO: Assignment Task 7: write function body
        self.mouse_pressed_pos = [event.xdata, event.ydata]
        self.closest_k = self.charges.get_closest(self.mouse_pressed_pos, 0.2)
        # End of Task 7; proceed to task 8.

    def on_mouse_release(self, event):
        ''' perform required actions when the mouse button is released

        If a charge should be deleted, delete the charge.
        If a charge should be added, add the charge.

        In all cases, redraw the figure with 32 fieldlines per charge
        and reset attributes as appropriate.
        '''
        # TODO: Assignment Task 8: write function body
        # resetting the values
        if self.closest_k is not None:
            # since there is a closest charge, it either needs to be moved
            # or deleted, depending on whether it is being dragged or not
            if self.dragging:
                self.charges.set_position(self.closest_k, [event.xdata, event.ydata])
            else:
                self.charges.delete_charge(self.closest_k)
        else:
            # not dragging, and didn't click anything, therefore add a new one
            if self.qadd != 0: # making sure you can't add neutral charges by accident
                self.charges.add_charge(self.qadd, self.mouse_pressed_pos)
                # reset the values of the attributes
                self.qadd = 0
                self.dragging = False
                self.closest_k = None
                self.mouse_pressed_pos = None
                # replot the charges and their field lines
                self.plot_fieldlines(self.field_lines_args[0], self.field_lines_args[1], self.field_lines_args[2], self.field_lines_args[3])
        # End of Task 8; no further tasks


app = None


# pylint: disable=C0111
def main():
    global app
    app = QApplication(sys.argv)
    w = MyMainWindow()
    w.show()
    app.exec()


if __name__ == '__main__':
    main()
