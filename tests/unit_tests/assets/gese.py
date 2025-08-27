from numpy import array, mean, diag


def forward_GeSe(pos, axes, **kwargs):
    Ge1, Ge2, Se1, Se2 = tuple(pos)
    a = axes[0, 0]
    b = axes[1, 1]
    x = mean([Ge1[0], Ge2[0] - 0.5])
    z1 = mean([Ge1[2], 1 - Ge2[2]])
    z2 = mean([Se2[2], 1 - Se1[2]])
    return [a, b, x, z1, z2]
# end def


def backward_GeSe(params, **kwargs):
    a, b, x, z1, z2 = tuple(params)
    Ge1 = [x, 0.25, z1]
    Ge2 = [x + 0.5, 0.75, 1 - z1]
    Se1 = [0.5, 0.25, 1 - z2]
    Se2 = [0.0, 0.75, z2]
    axes = diag([a, b, 20.0])
    pos = array([Ge1, Ge2, Se1, Se2])
    return pos, axes
# end def


# test GeSe monolayer
#                    a     b     x       z1       z2
params_GeSe = array([4.26, 3.95, 0.4140, 0.55600, 0.56000])
axes_GeSe_in = array([[4.26],
                      [3.95],
                      [20.0]])
pos_GeSe, axes_GeSe = backward_GeSe(params_GeSe)
elem_GeSe = 'Ge Ge Se Se'.split()


# random guess for testing purposes
hessian_GeSe = array([[1.00, 0.50, 40.0, 50.00, 60.0],
                      [0.50, 2.50, 20.0, 30.00, 10.0],
                      [40.0, 20.0, 70.0, 30.00, 10.0],
                      [50.0, 30.0, 30.0, 130.0, 90.0],
                      [60.0, 10.0, 10.0, 90.00, 210.0]])
