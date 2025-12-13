::: center
soilice: A modelling framework to test alternative process
representations in a frozen soil hydrological model\
Andrew Ireson\
2025-12-12
:::

# Developing Governing Equations for soilice

## Internal energy

Consider a control volume of water, of mass $M$ (kg) and temperature $T$
($^\circ$C). We define internal energy associated with sensible heat,
that is the sensible internal energy, $U_s$ (J), as

$$
U_s=M c_pT$$

where $c_p$ (J kg$^{-1}$ K$^{-1}$) is the specific heat capacity (which
is defined by this equation). The differential equation for the change
in sensible internal energy with time is

$$\frac{dU_s}{dt}=\frac{d}{dt}(M c_p T)=M c_p \frac{dT}{dt} + c_p T \frac{dM}{dt}$$

Latent heat is a form of internal energy, called latent internal energy,
$U_l$ (J), associated with the chemical bonds between molecules in a
particular phase (solid, liquid, gas) of a substance. We assume latent
heat is negative in the solid phase, zero in the liquid phase and
positive in the vapor phase, hence

$$\begin{array}{c}
    U_{li}=-M_i\lambda_f
    \\
    U_{ll}=0
    \\
    U_{lv}=-M_v\lambda_v
    \end{array}$$

where $\lambda_f$, the amount of heat needed to convert 1 kg of ice to 1
kg of liquid water, is a constant $\approx 3.34\times10^5$ J kg$^{-1}$
and $\lambda_v=2.50\times10^6$ J kg$^{-1}$ at $T=0^\circ$C).

Using this convention, the total latent internal energy is

$$
U_l=-M_i \lambda_f +M_v \lambda_v$$

and the change in latent internal energy is

$$\frac{dU_l}{dt}=-\lambda_f\frac{dM_i}{dt}+\lambda_v\frac{dM_v}{dt}$$

The total internal energy that concerns us in soils is given by the sum
of sensible and latent heat, i.e.

$$U=U_s+U_l$$

In a control volume containing ice, liquid water, vapor and soil solids
we have

$$U=U_{si}+U_{sl}+U_{sv}+U_{ss}+U_{li}+U_{ll}+U_{lv}$$

Adopting our conventions above we can drop $U_{ll}$.

## Mass and energy balance equations

We express the mass and internal energy per unit volume, i.e., as
$\bar{m}$ (kg m$^{-3}$) and $\bar{u}$ (J m$^{-3}$). The mass and energy
balance equations are

$$
\frac{d\bar{m}}{dt}=\frac{d\bar{m}_l}{dt}+\frac{d\bar{m}_i}{dt}=-\nabla q_l\rho_l$$

and

$$\frac{d\bar{u}}{dt}=
\frac{d\bar{u_{si}}}{dt}+
\frac{d\bar{u_{sl}}}{dt}+
\frac{d\bar{u_{ss}}}{dt}+
\frac{d\bar{u_{li}}}{dt}=
-\nabla (j_a+j_c+j_o)$$

where $\bar{u}_{ss}$ (J m$^{-3}$) is the sensible internal energy of the
soil solids, $j_a$, $j_c$ and $j_o$ (J m$^{-2}$ s$^{-1}$) are the
advective, conductive, and other heat flux components.

Expanding the left hand side of the energy balance equation gives

$$
    \begin{array}{ll}
        \displaystyle
\frac{d\bar{u}}{dt}
&
        \displaystyle
=\left(\bar{m}_{i} c_{pi} + \bar{m}_{l} c_{pl} + \bar{m}_{s} c_{ps}\right)
\frac{dT}{dt}
+c_{pl}T
\frac{d\bar{m}_{l}}{dt}
+(c_{pi}T-\lambda_f)
\frac{d\bar{m}_i}{dt}\\
    \end{array}$$

We define the bulk heat capacity, $c_{pb}$ (J kg$^{-1}$ K$^{-1}$) as

$$c_{pb}=\bar{m}_{i} c_{pi} + \bar{m}_{l} c_{pl} + \bar{m}_{s} c_{ps}$$

so we have

$$
    \begin{array}{ll}
        \displaystyle
\frac{d\bar{u}}{dt}
&
        \displaystyle
=c_{pb}
\frac{dT}{dt}
+c_{pl}T
\frac{d\bar{m}_{l}}{dt}
+(c_{pi}T-\lambda_f)
\frac{d\bar{m}_i}{dt}\\
    \end{array}$$

Now we want to eliminate the terms $d\bar{m}_l/dt$ and $d\bar{m}_i/dt$.
Under equilibrium conditions (which is a big assumption) a soil control
volume should reach a unique state on the basis of the temperature and
the total water. We ignore any hysteresis effects. We define the total
water content, $\theta_T$ (-) as being given by

$$\theta_T=\theta_L+\frac{\rho_I}{\rho_L}\theta_I$$

We defined $\psi_e$ (m) as the effective matric potential, which is
associated with $\theta_T$ by the soil characteristic curve
relationship,

$$\theta_T=M(\psi_e)$$

$\psi_f$ (m) is the matric potential that corresponds to the freezing
threshold, as defined by the GCE, where (Amankwah et al., 2021):

$$\psi_f=\frac{\lambda_f}{g}\ln\left(\frac{T+T_0}{T_0}\right)$$

where $T$ is the soil temperature in $^\circ$C, and $T_0 = 273.15$. If
$\psi_e\leq\psi_f$ there is no ice in the soil pore space. Now we say

$$\begin{array}{llll}
        \theta_T=M(\psi_e); & \theta_L=\theta_T; & \theta_I = 0; & \mathrm{if}\hspace{6pt} \psi_e\leq\psi_f\\
        \\
        \theta_T=M(\psi_e); & \theta_L=M(\psi_f(T))=F(T); & \theta_I = \frac{\rho_L}{\rho_I}(\theta_T-\theta_L); & \mathrm{if}\hspace{6pt} \psi_e>\psi_f\\
    \end{array}$$

The function $F(T)$ defines the freezing characteristic curve, and is
simply the combination of the GCE and the $\theta=M(\psi)$ relationship.
$M(\psi)$ and $F(T)$ can be differentiated, giving

$$\frac{d\theta_T}{d\psi_e}=M'(\psi_e)$$

and

$$\frac{d\theta_{L}}{dT}=G'(T)$$

So we can further say

$$
    \begin{array}{llll}
        \displaystyle \frac{d\theta_T}{dt} = M'(\psi_e)\frac{d\psi_e}{dt};
        & \displaystyle \frac{d\theta_L}{dt} = \frac{d\theta_T}{dt};
        & \displaystyle \frac{d\theta_I}{dt} = 0;
        & \displaystyle \mathrm{if}\hspace{6pt} \psi_e\leq\psi_f\\
        \\
        \displaystyle \frac{d\theta_T}{dt} = M'(\psi_e)\frac{d\psi_e}{dt};
        & \displaystyle \frac{d\theta_L}{dt} = F'(T)\frac{dT}{dt};
        & \displaystyle \frac{d\theta_I}{dt} = \frac{\rho_L}{\rho_I}\left(\frac{d\theta_T}{dt}-\frac{d\theta_L}{dt}\right);
        & \displaystyle \mathrm{if}\hspace{6pt} \psi_e>\psi_f\\
    \end{array}$$

Now, consider first the case where $\psi_e \leq \psi_f$ (i.e. the soil
is unfrozen). We have

$$\begin{array}{cc}
\displaystyle \frac{d\bar{m}_l}{dt}=\frac{d\bar{m}}{dt};
& 
\displaystyle \frac{d\bar{m}_i}{dt}=0
    \end{array}$$ 
    
Substituting this into Equation we have

$$
\frac{dT}{dt}
=
\frac{\displaystyle\frac{d\bar{u}}{dt}
-c_{pl}T
\frac{d\bar{m}}{dt}
}{c_{pb}
}$$

Next consider the case where $\psi_e > \psi_f$ (i.e. the soil is
frozen). This time we have

$$\frac{d\bar{m}_l}{dt}=\rho_lF'(T)\frac{dT}{dt}$$

and

$$\frac{d\bar{m}_i}{dt}=\frac{d\bar{m}}{dt}-\frac{d\bar{m}_l}{dt}
=
\frac{d\bar{m}}{dt}-\rho_lF'(\bar{m},T)\frac{dT}{dt}$$

Substituting these expressions into the Equation
[\[eq: UBal3DP\]](#eq: UBal3DP){reference-type="ref"
reference="eq: UBal3DP"} we get

$$\frac{d\bar{u}}{dt}
=
(c_{pb}
+\rho_l F'(T)
((c_{pl}-c_{pi})T+\lambda_f)
)\frac{dT}{dt}
+(c_{pi}T-\lambda_f)
\frac{d\bar{m}}{dt}$$

Rearranging we have

$$
\frac{dT}{dt}
=
\frac{\displaystyle
\frac{d\bar{u}}{dt}
-(c_{pi}T-\lambda_f)
\frac{d\bar{m}}{dt}
}{c_{pb}
+\rho_l F'(T)
((c_{pl}-c_{pi})T+\lambda_f)
}$$

which may then be solved by replacing the derivative terms on the right
hand side with the net fluxes.

## Finite difference formulation

For a one-dimensional vertical model, with $z$ (m) representing the
depth below the ground surface, it is good practice to use a
block-centred grid, where $z_n$ denotes the midpoint of a cell, and
$z_{n-1/2}$ and $z_{n+1/2}$ are the depths of the top and bottom faces
of the cell, respectively. Considering a single control volume at depth
index $n$, we can write

$$
    \begin{array}{rcl}
        \displaystyle
\left . \frac{d\bar{m}_l}{dt}\right |_n
+ \left . \frac{d\bar{m}_i}{dt}\right |_n
& 
\displaystyle 
= -\rho_l\nabla q_n
%&
%       \displaystyle
%=\frac{(q_{l,n-1/2}-q_{l,n+1/2})\rho_l}{z_{n+1/2}-z_{n-1/2}}
\\
        \displaystyle
\left . \frac{d\bar{u}}{dt}\right |_n
& 
\displaystyle 
= -\nabla j_n
%&
%       \displaystyle
%=\frac{(j_{a,n-1/2}-j_{a,n+1/2})
%+(j_{c,n-1/2}-j_{c,n+1/2})
%+(j_{o,n-1/2}-j_{o,n+1/2})
%}{z_{n+1/2}-z_{n-1/2}}
    \end{array}$$

where $\nabla q_n$ and $\nabla j_n$ are the net fluxes of mass and
energy into the control volume. Expressing these equations with total
mass, $\bar{m}$ and temperature $T$ as the dependent variables we have

$$\left . \frac{d\bar{m}}{dt} \right |_n
=-\rho_l \nabla q_n$$

$$\left . \frac{dT}{dt} \right |_n
=
\left \{
    \begin{array}{ll}
\displaystyle
\frac{ \nabla j_n
-c_{pl}T_n
\left . \frac{d\bar{m}}{dt} \right |_n
}{c_{pb}
}
,
&\displaystyle  \psi_{e,n} \leq \psi_{f,n}
\\
\displaystyle
\frac{\nabla j_n
-(c_{pi}T_n-\lambda_f)
\left . \frac{d\bar{m}}{dt} \right |_n
}{c_{pb}
+\rho_l f'(\bar{m_n},T_n)
((c_{pl}-c_{pi})T_n+\lambda_f)
}
,
&\displaystyle \psi_{e,n} > \psi_{f,n}
    \end{array}
\right .$$

Expanding the flux terms we have

$$
    \begin{array}{rcl}
\displaystyle 
-\rho_l\nabla q_n
&
        \displaystyle
=\frac{(q_{l,n-1/2}-q_{l,n+1/2})\rho_l}{z_{n+1/2}-z_{n-1/2}}
\\
\displaystyle 
-\nabla j_n
&
        \displaystyle
=\frac{(j_{a,n-1/2}-j_{a,n+1/2})
+(j_{c,n-1/2}-j_{c,n+1/2})
+(j_{o,n-1/2}-j_{o,n+1/2})
}{z_{n+1/2}-z_{n-1/2}}
    \end{array}$$

Here the individual flux terms are given by

$$\begin{array}{rl}
    \displaystyle
q_{l,n-1/2}&=-K_{n-1/2}\left(\frac{\psi_n-\psi_{n-1}}{z_i-z_{n-1}}-1\right)
\\
    \displaystyle
q_{l,n+1/2}&=-K_{n+1/2}\left(\frac{\psi_{n+1}-\psi_{n}}{z_{n+1}-z_{n}}-1\right)
\\
    \displaystyle
j_{a,n-1/2}&=q_{l,n-1/2}\rho_l c_{pl} T_{n-1}
\\
    \displaystyle
j_{a,n+1/2}&=q_{l,n+1/2}\rho_l c_{pl} T_{n}
\\
    \displaystyle
j_{c,n-1/2}&=-\kappa_{n-1/2}\left(\frac{T_n-T_{n-1}}{z_n-z_{n-1}}\right)
\\
    \displaystyle
j_{c,n+1/2}&=-\kappa_{n+1/2}\left(\frac{T_{n+1}-T_{n}}{z_{n+1}-z_{n}}\right)
\end{array}$$

For the advected heat fluxes here we assume the water flux is positive
in the $z$ direction, such that $T_{n-1}$ is the downstream temperature
at point $n-1/2$, and $T_n$ is the upstream temperature at point
$n+1/2$. Alternative formulations for advection may be considered, but
we will not expand on this here. The term $j_o$ will be zero for all $n$
except $n=1$, that is the upper boundary condition (i.e.
$j_{o,1/2}\neq 0$).

Consider a soil profile subject to a vertical infiltration flux $q_a$ (m
s$^{-1}$), a drainage flux, $q_d$ (m s$^{-1}$) and an evaporative loss
term $q_e$ (m s$^{-1}$). The following boundary conditions can be used,
where index $n=1$ represents the ground surface and $n=N$ is the lower
boundary:

$$
    \begin{array}{rl}
q_{l,1/2}& =q_a-q_e
\\
q_{l,N+1/2}& =q_d
\\
j_{a,1/2}& =j_{a,T}
\\
j_{c,1/2}& =j_{g,T}
\\
j_{a,N+1/2}& =j_{a,B}
\\
j_{c,N+1/2}& =j_{g,B}
\\
j_{o,1/2}& =j_r-j_h-j_e
\\
j_{o,N+1/2}& =0
    \end{array}$$

where all $j$ terms have units (J m$^{-2}$ s$^{-1}$), and $j_{a,T}$ is
advection with net infiltrating water, $j_{g,T}$ is the ground heat
flux, driven by conduction across the soil surface (important when there
is snow on the ground), $j_r$ is net radiation, $j_h$ is the sensible
heat lost (advection with the movement of air), $j_e$ is the latent heat
lost to evaporation (advection of water vapor), $j_{a,B}$ is advection
with water draining out the base of the soil and $j_{g,B}$ is the
conductive heat flux across the lower boundary. The different flux terms
in Eq [\[eq: BCs\]](#eq: BCs){reference-type="ref" reference="eq: BCs"}
may be zero under different scenarios (e.g. snow cover, bare soil, etc)
or assumptions.

We have here a system of ordinary differential equations which we solve
using an ODE solver with the method of lines.

## Constituitive relationships

Now, we need expressions for $M$, $M'$, $F$ and $F'$, all of which we
can get from the van Genuchten equations combined with the GCE equation.

$$S_e=(1+(\alpha \psi)^n)^{-m}$$

$$\theta=M(\psi)=\theta_r+(\theta_s-\theta_r)S_e$$

$$\frac{d\theta}{d\psi}=M'(\psi)=
\frac{-\alpha m (\theta_s-\theta_r)}{1-m}
S_e^{1/m}\left(1-S_e^{1/m}\right)^m$$

The GCE:

$$\psi_f(T)=\frac{\lambda_f}{g}\ln\left(\frac{T+T_0}{T_0}\right)\approx\frac{L_f}{g}\left(\frac{T}{T_0}\right)$$

and

$$\begin{array}{ll}
    F(T) & =M\left(\psi_f(T)\right)
    % = \theta_r+(\theta_s-\theta_r)\left(1+(\alpha\psi_f)^n\right)^{-m}
    \\
    & =
    \theta_r+(\theta_s-\theta_r)\left(1+\left(
    \frac{\alpha L_fT}{gT_0}
    \right)^n\right)^{-m}
\end{array}$$

For the derivative of the SFC let us define

$$S_f=\left(1+\left(\frac{\alpha L_f T}{gT_0}\right)^n
\right)^{-m}$$

Then

$$F'(T)=
\frac{-\alpha L_f m (\theta_s-\theta_r)}{gT_0(1-m)}
S_f^{1/m}\left(1-S_f^{1/m}\right)^m$$

Finally we need to define the hydraulic conductivity. In unfrozen
conditions this is

$$K(\psi)=K_sS_e^{1/2}\left(1-\left(1-S_e^{1/m}\right)^m\right)$$

In unfrozen conditions, we are still seeking an ideal relationship for
$K(\psi_e,T)$. For now we use impedence.
