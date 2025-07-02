
def create_connecting_line_seg(pt1, pt2):
      Plt_line_seg = parametric_plot3d((lambda lam: lam*pt1[0]+(1-lam)*pt2[0], lambda lam: lam*pt1[1]+(1-lam)*pt2[1], lambda lam: lam*pt1[2]+(1-lam)*pt2[2]), (0,1), frame=False, color='green', thickness=5)
      return Plt_line_seg



def f_param_fct_1(t):
     return t^2
def f_param_fct_2(t):
        return( t^3)
   
    
def f_param_fct_3(t):
     return (t^4)

param_fct_f = (f_param_fct_1, f_param_fct_2, f_param_fct_3)

Plot1 = parametric_plot3d(param_fct_f, (-sqrt(2.),sqrt(2.)), frame=False, color='red', thickness=5)

show(Plot1)

x1, x2, x3 = var( 'x1', 'x2', 'x3')


EquationXv = 4*x1^3*x2^2 + 27*x1*x2^4 - 16*x1^4*x3-144*x1*x2^2*x3 + 128*x1^2*x3^2 - 256*x3^2 == 0

Plot2 = implicit_plot3d(EquationXv, (x1,0,5), (x2,-2,5), (x3,-2,7), frame=False, color='blue', opacity=0.9)

show(Plot1 + Plot2)



t_val = 0.9
pt_X = vector([t_val^2, t_val^3, t_val^4])
pt_X_v = vector([t_val^2*4, t_val^3*(4), t_val^4*(-1)])
crit_pt = pt_X + pt_X_v

c = 0.2
v = (-4*c, 4*c, -1*c)
v = (-3*c, -2*c, 0*c)

Plot3 = parametric_plot3d((lambda lam: t_val^2*(1+lam*v[0]), lambda lam: t_val^3*(1+lam*v[1]), lambda lam: t_val^4*(1+lam*v[2])), (0,1), frame=False, color='green', thickness=5)
Plot4 = parametric_plot3d((lambda lam: t_val^2*(lam*1+v[0]), lambda lam: t_val^3*(lam*1+v[1]), lambda lam: t_val^4*(lam*1+v[2])), (0,1), frame=False, color='green', thickness=5)


show(Plot1 + Plot2 + Plot3+Plot4)








PlotLin1 = parametric_plot3d((lambda lam: lam*2, lambda lam: lam*3, lambda lam: lam*4), (-0.1,0.3), frame=False, color='red', thickness=5)


transl_plot_Lin2 = (-2,3,1)
Ray_vec1 = (1, 0, 0)
Ray_vec2 = (0, 1, 0)
Ray_vec3 = (0, 0, 1)
Ray_vec4 = (-1, -1, -1)


PlotVec1 = parametric_plot3d((lambda lam, mu: lam*2+mu*Ray_vec1[0]+transl_plot_Lin2[0], lambda lam, mu: lam*3+mu*Ray_vec1[1]+transl_plot_Lin2[1], lambda lam, mu: lam*4+mu*Ray_vec1[2]+transl_plot_Lin2[2]), (-0.3,0.3), (0,1), frame=False, color='blue', opacity=0.9)
PlotVec2 = parametric_plot3d((lambda lam, mu: lam*2+mu*Ray_vec2[0]+transl_plot_Lin2[0], lambda lam, mu: lam*3+mu*Ray_vec2[1]+transl_plot_Lin2[1], lambda lam, mu: lam*4+mu*Ray_vec2[2]+transl_plot_Lin2[2]), (-0.3,0.3), (0,1), frame=False, color='blue', opacity=0.9)
PlotVec3 = parametric_plot3d((lambda lam, mu: lam*2+mu*Ray_vec3[0]+transl_plot_Lin2[0], lambda lam, mu: lam*3+mu*Ray_vec3[1]+transl_plot_Lin2[1], lambda lam, mu: lam*4+mu*Ray_vec3[2]+transl_plot_Lin2[2]), (-0.3,0.3), (0,1), frame=False, color='blue', opacity=0.9)

x1 = (0.15*2+0.4*Ray_vec3[0]+transl_plot_Lin2[0], 0.15*3+0.4*Ray_vec3[1]+transl_plot_Lin2[1], 0.15*4+0.4*Ray_vec3[2]+transl_plot_Lin2[2])
x2 = (0.15*2, 0.15*3, 0.15*4)
u = ((x1[0]+x2[0]-1.7)/2, (x1[1]+x2[1]-0.4)/2, (x1[2]+x2[2]+1.7)/2)

LineSeg1 = create_connecting_line_seg(x1, u)
LineSeg2 = create_connecting_line_seg(x2, u)

show(PlotLin1+PlotVec1+PlotVec2+PlotVec3+LineSeg1+LineSeg2)