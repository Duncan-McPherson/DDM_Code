function domga_dt = Drive_ODE(t,omga,u, J, B, Tp_stat, Om1p, Tp_coul,   Om2p,   Tp_visc,   Tn_stat, Om1n, Tn_coul,   Om2n,   Tn_visc)

if abs(omga)<0.00001
    if (u==0)  Tf=0;  end
    if (u>0 && u<Tp_stat)  Tf=u;  end
    if (u<0 && u>Tn_stat)  Tf=u;  end
    if (u>0 && u>Tp_stat)  Tf=Tp_stat;  end
    if (u<0 && u<Tn_stat)  Tf=Tn_stat;  end
elseif omga>0
        Tf=Tp_stat*exp(-omga/Om1p)+Tp_coul*(1-exp(-omga/Om2p))+Tp_visc*omga;
elseif omga<0
    Tf=Tn_stat*exp(-omga/Om1n)+Tn_coul*(1-exp(-omga/Om2n))+Tn_visc*omga;
end

domga_dt=-B*omga/J+   u/J     - Tf/J;