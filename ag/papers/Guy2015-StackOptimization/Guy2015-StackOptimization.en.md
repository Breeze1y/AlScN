<!-- PDF page 1; printed page 3482 -->

3482 IEEE TRANSACTIONS ON ELECTRON DEVICES, VOL. 62, NO. 11, NOVEMBER 2015
Investigation of Forming, SET, and Data Retention
of Conductive-Bridge Random-Access
Memory for Stack Optimization
Jérémy Guy, Gabriel Molas, Philippe Blaise, Mathieu Bernard, Anne Roule, Gilles Le Carval, Vincent Delaye,
Alain Toffoli, Gérard Ghibaudo, Fellow, IEEE, Fabien Clermidy, Barbara De Salvo, and Luca Perniola
Abstract— In this paper, we investigate in depth Forming,
SET, and Retention of conductive-bridge random-access
memory (CBRAM). A kinetic Monte Carlo model of the CBRAM
has been developed considering ionic hopping and chemical
reaction dynamics. Based on inputs from ab initio calculations
and the physical properties of the materials, the model offers the
simulation of both the Forming/SET and the Data Retention
operations. It aims to create a bond between the physics at
atomic level and the device behavior. From the model and
experimental results obtained on decananometric devices, we
propose an understanding of the physical mechanisms involved
in the CBRAM operations. Using the consistent Forming/SET
and Data Retention model, we obtained good agreement with
the experimental data. Finally, the impact of each layer of the
CBRAM on the Forming/SET behavior is decorrelated, allowing
an optimization of the performance.
Index Terms — Conductive-bridge random-access memory
(CBRAM), Data Retention, Forming, kinetic Monte Carlo (KMC)
modeling, resistive memory, resistive RAM, SET, transition-state
theory (TST).
I. I NTRODUCTION
C
ONDUCTIVE-BRIDGE random-access memory
(CBRAM) is a promising nonvolatile memory
technology envisaged as an alternative to ﬂash memory. The
strength of this technology resides in its low power
consumption, high speed, simple structure, and ease of
integration in the back end of a logic process. It has
been accepted that this nonvolatile memory is based on
the reversible formation and dissolution of a conductive
Manuscript received March 6, 2015; revised August 13, 2015; accepted
September 1, 2015. Date of publicatio n September 25, 2015; date of current
version October 20, 2015. The review of this paper was arranged by Editor
J. S. Suehle. (Corresponding author: Jérémy Guy.)
J. Guy, G. Molas, P. Blaise, M. Bernard, A. Roule, G. Le Carval,
V . Delaye, A. Toffoli, F. Clermidy, B. De Salvo, and L. Perniola are with
Commisariat à l’Énergie Atomique –Laboratoire d’Électronique et de
Technologie de l’Information, Gre noble 38054, France (e-mail: jeremy.
guy@cea.fr; gabriel.molas@cea.fr; philippe.blaise@cea.fr; mathieu.bernard
@cea.fr; anne.roule@cea.fr; gilles.lecar val@cea.fr; vicent.delaye@cea.fr;
alain.toffoli@cea.fr; fabien.clermidy@cea.fr; barbara.desalvo@cea.frl;
luca.perniola@cea.fr).
G. Ghibaudo is with Institut de Microélectronique, Électromag-
nétisme et Photonique–Laboratoire d’Hyperfréquences et de Caractérisation,
Minatec/Grenoble Institute of Tec hnology, Grenoble 38016, France (e-mail:
ghibaudo@minatec.inpg.fr).
Color versions of one or more of the ﬁgures in this paper are available
online at http://ieeexplore.ieee.org.
Digital Object Identiﬁer 10.1109/TED.2015.2476825
ﬁlament (CF) inside a resistive layer (RL) [1], [2]. Various
RL materials have been studied from chalcogenides (GeSe and
so on [3], [4]) to oxides (Al
2O3 and so on [5], [6]), these latter
providing an improved Data Retention behavior. However, the
main drawback of oxide-based CBRAM is the necessity
of a high-voltage Forming operation. To respond to the
electrical needs of this unique Forming cycle, a complex
circuitry may be required at the cost of an integration density
loss [7]. Thus, the technical knowledge of the CBRAM keeps
growing, leading to high-density integration [8], whereas
the physical understanding of the mechanisms involved
during the SET, RESET, and Retention still has to be
clearly understood. Several studies have been carried out in
order to investigate these mechanisms, either by focusing
on CF as the growing entity at the device scale, based on
continuous [9] or discrete phenomena [10], or by performing
the ﬁrst principle calculations at the materials level to get
insights into the ﬁlament composition [11]. In this paper, we
present and describe an extended version of the kinetic
Monte Carlo (KMC) simulation introduced in [12].
We propose a consistent model of Forming/SET and Data
Retention using a trap-assisted tunneling (TAT) electrical
conduction model. This KMC simulation creates a bridge
between the physics at the atomic scale and the device
behavior and performance. From ab initio inputs and material
physical properties, based on the atomic and electronic
migration, we propose a physical interpretation of the
mechanisms involved in the switching operation and
Retention. Moreover, an experimental study on nanotrench
CBRAM structures [6], with various material combinations,
has been conducted and coupled to the KMC simulations
to extract the impact of each component on the CBRAM
characteristics. A speciﬁc attention is given to the time–
voltage dilemma, investigating the time–voltage dependence
of the Forming and SET characteristics. Finally, the material
properties and the understanding of the mechanisms involved
in the CBRAM operation allowed us to put forward the
optimized CBRAM stacks for improved switching speeds.
II. T
ECHNOLOGICAL DETAILS
The nanotrench CBRAM structure studied in this paper
is represented (Fig. 1). The bottom electrode (BE) is
deﬁned by a metallic liner deposited by Chemical Vapor
0018-9383 © 2015 I EEE. Personal u se is perm itted, but republication/redistribution requires IEEE permission.
See http://www.ieee.org/publications_standards/publications/rights/index.html for more information.
Authorized licensed use limited to: BEIJING INSTITUTE OF TECHNOLOGY. Downloaded on January 05,2026 at 01:22:27 UTC from IEEE Xplore.  Restrictions apply. 

---

<!-- PDF page 2; printed page 3483 -->

GUY et al.: INVESTIGATION OF FORMING, SET, AND DATA RETENTION OF CBRAM 3483
Fig. 1. Schematic process ﬂow description of the nanotrench CBRAM and
transmission electron microscopy cross sections.
TABLE I
STUDIED SAMPLES AND MATERIAL STACK VARIATIONS
Deposition in an etched via. The thickness of the layer,
down to 5 nm, controls the length of the BE. The via
is then ﬁlled with SiO 2, and the BE contact is opened
by chemical mechanical plan arization to create a ring.
The second critical dimension of the memory is deﬁned
by a nanotrench, etched in a SiN capping layer deposited
over the planarized BE. This nanotrench is deﬁned by
electron-beam lithography down to 50 nm and leads to the
reduction of the CBRAM active area down to 5 × 50 nm
2.
A metal oxide acting as an RL (Al 2O3 deposited by ALD
or a metal oxide MO x deposited by Physical Vapor Deposi-
tion (PVD)) is then deposited in the nanotrench and capped by
a PVD Cu-based top electrode (TE). Table I shows the samples
(S1–S6) with materials stack variations, studied to understand
the impact of each layer on the CBRAM operation.
III. M
ODEL DESCRIPTION
A. Device Operation Principle
The CBRAM tec hnology is based on the reversible
formation and dissolution of a metallic CF inducing a change
in the resistance state of the cell. This CF has been exper-
imentally observed with conductive-AFM tomography [13].
The formation and dissolution process is known to be related
to the oxidation and reduction of the metal composing the
TE or BE [1], [2]. In this paper, the Cu-based TE is considered
as the metal source for the CF creation. In order to initially
create the CF in the ﬁrst cycle, a high Forming voltage is
required. This step allows the metallic ions to ﬂow through the
RL toward the opposite electrode and to create the CF. Similar
to the SET operation and in opposition to the RESET, Forming
consists in switching the cell from a high resistive state (HRS)
to a low resistive state (LRS). After the Forming step, the
cell operates at lower voltages. Once the data are stored, the
main concern stands in its retention over time. Although the
CBRAM is considered as a nonvolatile memory, the CF tends
to dissolve itself over time, which ultimately leads to the loss
of data. The CF dissolution rate is activated in temperature.
Both the LRS and the HRS are potentially unstable states and
Fig. 2. Particle transition mechanisms with the activation energy modiﬁed
by electric ﬁeld and energy level difference. (a)–(d) E A is the energy barrier
height, q is the electrical charge of the hopping entity, V is the electrical
potential, and dhopping is the hopping distance.
Fig. 3. Simulated physical mechanisms considered in the KMC environment.
have to be carefully studied to evaluate the reliability of the
memory device [6].
B. KMC Environment
The aim of our simulation framework is to create a con-
sistent environment, in which Forming, SET, and Retention
can be simulated by sharing the same physical parameters
of the materials. The model is based on atomic migration
and oxidation reduction reactions treated at the atomic and
electronic level. Each and every atom is treated individually,
according to their electrical an d thermal environment. The
CBRAM operation results in a succession of various physical
mechanisms which ultimately leads to a change in resistance
state. To simulate the CBRAM operation, the KMC model
relies on a 2-D grid corresponding to every possible atomic
and ionic positions inside the RL. The grid parameter (hopping
distance between two adjacent points) is deﬁned by the atomic
Cu–Cu distance in the RL. Some insights into this latter can be
obtained by the ab initio calculations [12], [14] and depends
on the RL.
C. Physical Mechanisms and KMC Interpretations
1) Redox Reaction and Ionic Migration: Both redox
reactions and ionic migration can be correlated with a particle
transitioning from an energy state to a more favorable one,
overcoming an energy barrier [15] (Fig. 2). The oxidation
corresponds to the transfer of an electron from a Cu atom at
the edge of the RL to a neighboring metal atom. Similarly,
the reduction is considered as the transfer of an electron
to a Cu
+ ionized atom from the BE or from the growing
ﬁlament. Fig. 3 shows the different mechanisms previously
stated and the environment, in which they take place. The
Authorized licensed use limited to: BEIJING INSTITUTE OF TECHNOLOGY. Downloaded on January 05,2026 at 01:22:27 UTC from IEEE Xplore.  Restrictions apply. 

---

<!-- PDF page 3; printed page 3484 -->

3484 IEEE TRANSACTIONS ON ELECTRON DEVICES, VOL. 62, NO. 11, NOVEMBER 2015
particles hopping and redox reactions are described and
modeled through the transiti on-state theory (TST),
which introduces the mechanism rates [16], under the
name of Eyring–Polanyi equation, in agreement with the
Maxwell–Boltzmann statistic
Γ = v exp
(
− E
A
kB T
)
(1)
where Γ is the mechanism rate, v is the vibration constant,
EA is the energy barrier height, kB is the Boltzmann constant,
and T is the temperature. Considering the ﬁlament unformed,
the joule heating is neglected and the temperature is ﬁxed
by the external conditions. The mechanism rate expressed
in s−1 corresponds to the inverse of the required time for
a mechanism deﬁned by its activation energy EA to occur.
EA depends on the involved mechanism and can be modulated
(Fig. 2) by the applied electric ﬁeld as well as by the difference
of energy level between the two states, following :
E′
A = EA − Qd ∈
2 − Δ E
2 (2)
where E′
A is the modulated energy barrier height, Q is
the particle charge, d is the physical distance between the
two states, ∈ is the applied electric ﬁeld, and Δ E is the
energy difference between the two states. Moreover, it has
been demonstrated [17] that the local electric ﬁeld experienced
inside an oxide is modiﬁed by a polarization of the atomic
bonds. The local electric ﬁeld can be expressed as
∈
local =∈ applied (1 + Lχ) (3)
where L is the Lorentz factor, χ is the electric susceptibility
χ = εr − 1, and εr is the relative permittivity. By combin-
ing (1)–(3), it is possible to obtain the complete mechanism
rate equation for the CBRAM system
Γ = v exp
(
− E
A − 1
2(Qd(1 + Lχ) ∈− Δ E)
kB T
)
. (4)
In order to simplify, we ﬁxed v = 1013 s−1 [10] and
L = 1/3 [17]. We assume that the duration of each mechanism
is much shorter than the time it takes before it occurs. The
system can be considered as memoryless similar to atomic
exponential decay. This allows us to express the probability P
of each action for a time duration t as
P = 1 − exp(−Γ t). (5)
This equation bonds the physical parameters to the
occurring probability of each mechanism.
2) Electrical Conduction Model: As shown in (4), the
mechanisms are driven by the electric ﬁeld and the tempera-
ture. The temperature being ﬁxed by the external conditions,
we only have to calculate the e lectric ﬁeld on every points
of the grid. Moreover, the current ﬂowing through the cell
is another important matter, as it permits the monitoring of
the CF formation. In order to save computation time without
losing much accuracy, we chose to calculate at the same
time the current ﬂowing through the cell and the electrical
potential at every points of the grid. The current is calculated
Fig. 4. KMC ﬂowchart.
by solving the resistance grid between the atoms. We compute
the resistance between the two atoms using the TAT
R = R0 exp
(
2
√
 2m∗EC
ℏ (d − a)
)
(6)
where R0 is the quantum resistance (12.9 k /Omega1), m∗ is the
effective mass of electron for a speciﬁc RL, EC is the
energy barrier for electron tunneling, d is the distance between
the two atoms, and a is the contact distance between the
two atoms. The resistance depends on the distance between the
two atoms until they reach a contact position ( d = a) leading
to a quantum resistance. Using the atoms position, we create a
resistive grid between the atoms and solve the Kirchhoff law
to obtain the current ﬂowing through the cell and the voltage at
each atom location. Finally, with the ﬁnite difference method,
we acquire the voltage value across the RL. The TAT model
is in agreement with the previ ous statements [6] regarding
the progressive decrease in res istance while increasing the
Cu concentration into the RL.
D. KMC Operating
The KMC resolution follows the ﬂowchart in Fig. 4. It starts
by initializing the environment: grid, electrodes, and time.
Then, it calculates the current, voltage, and ﬁeld across the RL,
then it deduces the rates of the v arious possible mechanisms.
The Retention and Forming/SET regimes follow the same path,
with the only difference being the driving force. The main
driving force for the switching process is the voltage applied,
whereas none is applied during the Retention, and only the
temperature runs the process. Based on the computed mecha-
nisms rates, an optimized iteration duration t is determined
t =− ln(1 − rand)
ΣΓ (7)
where rand is a randomly generated value (ranging between
0 and 1) and corresponds to th e probability of at least one
event occurring, and ΣΓ is the sum of all the implied
mechanisms rates. From this iteration times, it computes all the
mechanisms probabilities following (5). After normalizing all
the probabilities, the occurring events are chosen depending
on their probability weight and the iteration time, then the
particles locations and states are updated. Finally, the current,
Authorized licensed use limited to: BEIJING INSTITUTE OF TECHNOLOGY. Downloaded on January 05,2026 at 01:22:27 UTC from IEEE Xplore.  Restrictions apply. 

---

<!-- PDF page 4; printed page 3485 -->

GUY et al.: INVESTIGATION OF FORMING, SET, AND DATA RETENTION OF CBRAM 3485
Fig. 5. Simulated ﬁlament growth and corresponding tForming extraction.
voltage, and ﬁeld across the RL are solved, and the model
determines if another loop is required in a function of the
switching or Retention conditions. In the end, the model is
able to produce a full simulate d data Retention or writing
process (Fig. 5).
IV . E
XPERIMENTAL RESULTS
A. Forming and SET Analyses
In this section, we analyze how each layer of the
CBRAM stack (TE, RL, and BE) impacts the Forming and
SET operations combining the ch aracterizations results and
the KMC simulation. The elect rical characterizations have
been performed using an Agilent parameter analyzer, offering
a good control of the limiting current (compliance current),
responsible for the LRS resistance. The Forming/SET voltages
were obtained using quasi-static programming conditions,
whereas the Forming/SET time were obtained using pulsed
programming and extracted as [6]. In regard to the KMC,
a SET can be simulated by creating a residual ﬁlament in
the RL and follows the same trend as a Forming with the
reduced initial RL thickness and resistance. A speciﬁc focus is
given to the relation between the Forming time and the voltage.
Indeed, the t
Forming (VPulse) characteristics show an exponential
dependence, meaning that increasing the programming pulse
drastically reduces the required Forming time. Moreover,
the slope of these characteristics illustrates the time–voltage
dilemma of resistive RAMs; a sharp slope indicating that short
programming times at high voltage and good disturb immunity
at low voltage may be combined.
1) Impact of the Top Electrode: Pure Cu as well as CuTe
x
TE, both with the same TiN BE were investigated, the
restrained Te concentration allows us to consider a unique
energy barrier height of Cu redox reactions [18], and ionic
hopping energy barrier height in the RL. The only variable we
consider by changing the TE is the metal work function (/Phi1
TE),
with, respectively, 5.2 and 5 eV for Cu and CuTe x .F i g .6
shows the tForming (VPulse) experimental and simulated curves.
We can see that an increase of /Phi1TE by 0.2 eV leads to a
tForming(VPulse) shift of 0.25 V toward lower voltages. The
TE work function has a double impact on the switching
operation (Fig. 6). Indeed, the difference of work function
between the TE and the BE induces a ﬂat band voltage VFB
Fig. 6. (a) Measured and simulated Forming time as a function of the applied
voltage for CBRAM integrating Cu or CuTe x (S1 versus S2—Table I) as the
TE (ion supply layer). (b) Schematics illustrating the impact of /Phi1TE on the
electric ﬁeld in the RL and on the redox reaction.
which is reduced as /Phi1TE increases
VFB = /Phi1BE − /Phi1TE. (8)
This ﬂat band voltage adds up to the applied voltage and
inﬂuences all the CBRAM physical mechanisms. Its impact
is related to the RL thickness. However, the Fermi level of
the TE corresponds to the ﬁnal state of the electrons during
oxidation. A higher /Phi1TE leads to a lower Fermi level and a
faster oxidation reaction in agreement with (4), independently
of the RL thickness. The two inﬂuences of TE are antagonist,
and the later one is dominant here as an increase of /Phi1
TE
induces a diminution of tForming.
2) Impact of the Bottom Electrode: For the BE impact,
we focused the study on two materials: TiN- and Si-doped W,
both with the same CuTe
x TE. These materials have two differ-
ent work functions, respectively, 4.8 and 4.4 eV . Fig. 7 shows
that an increase of /Phi1BE by 0.4 eV leads to a tForming(VPulse)
shift of 0.5 V toward lower voltages. Once again, the
electrode work function has two effects on the switching
operation (Fig. 7). First, a higher /Phi1BE leads to a lower
Fermi level (initial electron energy during reduction reaction)
and a slower reduction reaction of the ﬁrst ions reaching
the BE (4). Then, a higher /Phi1
BE increases VFB (8), which
accelerates the oxidation reacti on. Seeing that the increase
of /Phi1BE reduces tForming, we can conclude that the limiting
mechanisms in regard of the redox reaction are the oxidation
at the TE.
3) Impact of the Resistive Layer: T h eR Lh a sam u c hw i d e r
and complicated impact on th e switching process as many
mechanisms are involved. First, Fig. 8 shows an experimental
Authorized licensed use limited to: BEIJING INSTITUTE OF TECHNOLOGY. Downloaded on January 05,2026 at 01:22:27 UTC from IEEE Xplore.  Restrictions apply. 

---

<!-- PDF page 5; printed page 3486 -->

3486 IEEE TRANSACTIONS ON ELECTRON DEVICES, VOL. 62, NO. 11, NOVEMBER 2015
Fig. 7. (a) Measured and simulated Forming time as a function of the applied
voltage for CBRAM integrating WSi or TiN (S3 versus S2—Table I) as the
BE (inert electrode). (b) Schematics illustrating the impact of /Phi1BE on the
electric ﬁeld in the RL and on the redox reaction.
Fig. 8. (a) Measured and simulated Forming time as a function of the applied
voltage for the 3.5 nm (S2) and 5 nm (S4) RL thickness. (b) Schematics
illustrating the impact of the RL thickness on the electric ﬁeld.
shift and a change of slope of tForming(VPulse) toward lower
voltages by reducing the RL thickness, in agreement with
our simulations. An RL thickness reduction of 1.5 nm turns
into a reduction of Forming voltage between 1.25 V (for long
Fig. 9. (a) Measured and simulated Forming time as a function of the pulse
voltage for Al2O3-based (S4) and MO x -based (S5) CBRAM. (b) Schematics
illustrating the impact of /Phi1RL on the redox reaction and the hopping distance
difference.
Forming times of 100 ms) and 2 V (for short Forming times
of 100 ns) in the studied area. The electric ﬁeld is tied to
the voltage and distance by ∈= d−1 V . For the sake of
simpliﬁcation and to clarify this point, the distance d can be
assimilated to the RL thickness. This means that the electric
ﬁeld increases faster with voltage for a thinner RL. The
change of slope of t
Forming(VPulse) can thus be explained by an
enhancement of the voltage impact on the electric ﬁeld for a
thin RL. Moreover, the reduction of the RL thickness translates
into a shorter ionic migration path and a faster switching, in
part responsible for the curve shift. We have then compared
the Forming characteristics of Al
2O3 and MOx RL. A steeper
slope and a faster behavior are measured for MO x ,a ss h o w n
in Fig. 9. The redox reactions take place at the interface
between the RL and the electrodes. The RL plays a great role
in these reactions, and regarding the TST affects the initial and
ﬁnal energy levels of the electrons during the redox reactions.
Indeed, during the oxidation reaction, the electrons initial
energy state corresponds to the energy level of Cu atom traps
inside the RL. Moreover, regarding the reduction reaction, the
electrons ﬁnal energy state corresponds to the energy level of
Cu ions trap in the RL. The trap depths strongly depend on the
RL and can be extracted from the ab initio calculations [12].
In addition, the trap depths directly inﬂuence the reaction,
as the difference of energy levels between the initial and ﬁnal
electrons state accelerates the reaction (4). Thus, deep Cu atom
traps (initial state) induce a low energy difference between the
Authorized licensed use limited to: BEIJING INSTITUTE OF TECHNOLOGY. Downloaded on January 05,2026 at 01:22:27 UTC from IEEE Xplore.  Restrictions apply. 

---

<!-- PDF page 6; printed page 3487 -->

GUY et al.: INVESTIGATION OF FORMING, SET, AND DATA RETENTION OF CBRAM 3487
Fig. 10. Summary of the impact of each physical parameter of the CBRAM stack on the time–voltage Forming characteristics.
initial and ﬁnal electrons states and lead to a slow oxidation
reaction [Fig. 9(b)]. On the contrary, deep Cu ion traps
increase the energy states difference and increase the reduction
reaction speed. These trap depths can shift the t
Forming(VPulse)
curve but cannot explain the change of slope, as shown
in Fig. 9. To address the change of slope, we have to look
into parameters enhancing the electric ﬁeld. By referring to (4),
we see that the local electric ﬁel d is enhanced by the electric
permittivity as well as by the distance between two transition
states of the TST. Both of these parameters are different
for the two RL and participate to the change of slope: the
higher the permittivity and hoppi ng distance, the bigger the
slope. The distance between the two states can be extracted by
the ab initio calculation, as previously stated, and is affected
by the RL density. It corresponds to the average value of the
Cu–Cu distance in the RL for the ionic migration and in the
interfacial region for the redox reactions.
4) SET and Forming—Conclusion: It has been demon-
strated that by tuning the material stack and properties,
we are able to control the Forming and thus the SET operation
in two ways, as shown in Fig. 10. It is possible to shift the
t
Forming(VPulse) characteristic by tuning the energy levels of
the various layers (/Phi1RL, /Phi1TE,a n d/Phi1BE) but also to change the
characteristics slope by altering the RL features enhancing the
electric ﬁeld (thickness, perm ittivity, and hopping distance).
Both the shift and slope changes of the tForming(VPulse) char-
acteristics offer a way to incr ease the SET/Forming speed
(or to reduce the voltage), but in addition, the slope tuning
brings on a way to improve the disturb immunity. This disturb
immunity corresponds to the switching time at low voltages
(typically V
read) and can be tied to a better HRS Retention.
The disturb immunity is improved for high dielectric constant
Fig. 11. Calculation of the SET time/vo ltage dependence for various residual
ﬁlaments in the RL, corresponding to various initial ROFF values.
and thin RL. Finally, it is possible to extrapolate VSET from a
different initial ﬁlament height corresponding to various initial
HRS resistances, as shown in Fig. 11.
B. Data Retention Modeling
As previously stated, the goal of the KMC model is to create
a unique environment for Forming/SET and Data Retention.
This means that the observed behavior for the HRS and LRS
Retention is simulated with the same physical properties and
parameters as for the switching. The only difference resides
in the main accelerating factor, being the voltage for Forming/
SET and the temperature for Data Retention. Fig. 12 shows
the measured and simulated 100 °C Retention characteristics
for both the HRS and the LRS. The impact of the compliance
current used during the Forming/SET on the LRS Retention
capabilities can be seen in accordance with [6]: the higher the
compliance current, the lower the initial resistance, and the
Authorized licensed use limited to: BEIJING INSTITUTE OF TECHNOLOGY. Downloaded on January 05,2026 at 01:22:27 UTC from IEEE Xplore.  Restrictions apply. 

---

<!-- PDF page 7; printed page 3488 -->

3488 IEEE TRANSACTIONS ON ELECTRON DEVICES, VOL. 62, NO. 11, NOVEMBER 2015
Fig. 12. LRS and one HRS Data Retention experimental results and
KMC simulations showing the impact of the initial RON on the Retention
capabilities at 100 °C and the dissolution of the CF. The yellow curves
share the same LRS resistance. 1: initial HRS. 2: ﬁnal HRS. 3: ﬁnal LRS.
4: initial LRS. All corresponding to the data: R
init−2.
Fig. 13. LRS Data Retention experimental results and KMC simulations
showing the impact of the temperatu re in the Retention capabilities and the
dissolution of the CF.
better the Retention. This can be explained by the fact that
the core of the ﬁlament, responsible for the ohmic conduction,
remains intact longer for a thick ﬁlament, which enhances
its Retention capabilities. Then, HRS is described by the
presence of a residual ﬁlament, whose height depends on the
HRS resistance, in agreement w ith previous statements [6].
HRS resistance decreases over time as the dissolution of the
remaining ﬁlament tends to recreate a conducting path in
the RL. The reformed conductive path can be seen (Fig. 12)
on the ﬁnal HRS showing only little differences with the
ﬁnal LRS. Moreover, as stated before, the temperature pushes
the Retention operation. Fig. 13 shows the temperature impact
on the Retention with the corresponding ﬁnal states. We can
see an increase of resistance value with the time spent at
high temperature. The resistance increase corresponds to the
loss of data over time and is due to the CF dissolution. This
dissolution is much more visible for high temperature similar
to the resistance increase. The higher the temperature, the more
dissolved the CF, and the worst the Retention. It is possible to
extract an activation energy of 1.2 eV from the experimental
and simulated data. The KMC Retention simulation shows
great results in agreement with the studied devices, in term
of temperature and compliance current dependence on both
the LRS and the HRS (not shown here) [6]. Its ability to
predict the LRS and HRS behavior at various temperatures
offers a way to optimize the SET conditions to obtain the
best Retention window possible for the targeted application
and temperature. The programming routine of a CBRAM can
then be adapted to achieve either low-power programming and
moderate-temperature operation or higher programming power
and high-temperature operation [19].
V. C
ONCLUSION
We presented a KMC simulation able to recreate Forming,
SET operation, and Data Retention of a CBRAM at atomic
level. Good agreement with the experimental data obtained on
decanometric CBRAM devices has been shown. From the joint
study of modeling and electrical characterization, we described
the relationship between th e physical parameters of each
layer and the device performan ce. Therefore, we presented
a way to optimize the device in order to improve the time–
voltage dilemma. Finally, we illustrated the KMC simulation
ability to reproduce a Retention study from the simulated SET
depending on the temperature accelerating the process.
R
EFERENCES
[1] R. Waser and M. Aono, “Nanoionics-based resistive switching
memories,” Nature Mater., vol. 6, no. 11, pp. 833–840, Nov. 2007.
[2] R. Waser, R. Dittmann, G. Staikov, and K. Szot, “Redox-based
resistive switching memories—Nanoionic mechanisms, prospects, and
challenges,” Adv. Mater. , vol. 21, nos. 25–26, pp. 2632–2663,
Jul. 2009.
[3] R. Soni, M. Meier, A. Rüdiger, B. Holländer, C. Kügeler, and R. Waser,
“Integration of ‘Ge
xSe1−x’ in crossbar arrays for non-volatile memory
applications,” Microelectron. Eng., vol. 86, nos. 4–6, pp. 1054–1056,
Apr./Jun. 2009.
[4] J. Guy et al., “Impact of Sb doping on power consumption and retention
reliability of GeS 2 based conductive bridge random access memory,”
Thin Solid Films , vol. 563, pp. 15–19, Jul. 2014.
[5] S. Maikap, S. Z. Rahaman, T. Y . Wu, F. Chen, M. J. Kao, and M. J. Tsai,
“Low current (5 pA) resistive switching memory using high- k Ta2O5
solid electrolyte,” in Proc. ESSDERC, Sep. 2009, pp. 217–220.
[6] J. Guy et al. , “Investigation of the physical mechanisms governing
data-retention in down to 10 nm nano-trench Al 2O3/CuTeGe con-
ductive bridge RAM (CBRAM),” in IEDM Tech. Dig. , Dec. 2013,
pp. 30.2.1–30.2.4.
[7] S. Sills et al. , “A copper ReRAM cell for storage class memory
applications,” in VLSI Tech. Dig., Jun. 2014, pp. 1–2.
[8] J. Zahurak et al., “Process integration of a 27 nm, 16 Gb Cu ReRAM,”
in IEDM Tech. Dig., Dec. 2014, pp. 6.2.1–6.2.4.
[9] D. Ielmini, “Modeling the universal set/reset characteristics of bipolar
RRAM by ﬁeld- and temperature-driven ﬁlament growth,” IEEE Trans.
Electron Devices, vol. 58, no. 12, pp. 4309–4317, Dec. 2011.
[10] J. R. Jameson et al., “Effects of cooperative ionic motion on program-
ming kinetics of conductive-bridge memory cells,” Appl. Phys. Lett. ,
vol. 100, no. 2, p. 023505, 2012.
[11] K. Sankaran et al. , “Modeling of copper diffusion in amorphous alu-
minum oxide in CBRAM memory stack,” ECS Trans., vol. 45, no. 3,
pp. 317–330, 2012.
[12] J. Guy et al. , “Experimental and theore tical understanding of form-
ing, SET and RESET operations in conductive bridge ram (CBRAM)
for memory stack optimization,” in IEDM Tech. Dig. , Dec. 2014,
pp. 6.5.1–6.5.4.
[13] U. Celano, Y . Y . Chen, D. J. Wouters, G. Groeseneken, M. Jurczak,
and W. Vandervorst, “Filament observation in metal-oxide resistive
switching devices,” Appl. Phys. Lett. , vol. 102, no. 12, p. 121602,
2013.
[14] L. Goux et al. , “Field-driven ultrafast sub-ns programming in
WAl2O3TiCuTe-based 1T1R CBRAM system,” in VLSI Tech. Dig. ,
Jun. 2012, pp. 69–70.
[15] S. Larentis, F. Nardi, S. Balatti, D. C. Gilmer, and D. Ielmini, “Resistive
switching by voltage-driven ion mi gration in bipolar RRAM—Part II:
Modeling,” IEEE Trans. Electron Devices, vol. 59, no. 9, pp. 2468–2475,
Sep. 2012.
[16] K. J. Laidler and M. C. King, “Devel opment of transition-state theory,”
J. Phys. Chem. , vol. 87, no. 15, pp. 2657–2664, 1983.
Authorized licensed use limited to: BEIJING INSTITUTE OF TECHNOLOGY. Downloaded on January 05,2026 at 01:22:27 UTC from IEEE Xplore.  Restrictions apply. 

---

<!-- PDF page 8; printed page 3489 -->

GUY et al.: INVESTIGATION OF FORMING, SET, AND DATA RETENTION OF CBRAM 3489
[17] J. W. McPherson, R. B. Khamankar, and A. Shanware, “Complementary
model for intrinsic time-depende nt dielectric breakdown in SiO 2
dielectrics,” J. Appl. Phys. , vol. 88, no. 9, pp. 5351–5359,
Mar. 2000.
[18] F. De Stefano et al. , “Inﬂuence of metal electrode stoichiometry
on the electron barrier height at Cu xTe1−x/Al2O3 interfaces for
CBRAM applications,” Microelectron. Eng. , vol. 120, pp. 9–12,
May 2014.
[19] M. Barci et al., “Impact of SET and RESET conditions on CBRAM high
temperature data retention,” in Proc. Int. Rel. Phys. Symp. , Jun. 2014,
pp. 5E.3.1–5E.3.4.
Jérémy Guy received the Engineering degree in
material science from the Institut National des Sci-
ences Appliquées de Lyon (INSA), Lyon, France,
in 2012, and the M.S. degree in microelectronic and
embedded system from INSA, and the University
of Lyon, Lyon, in 2012. He is currently pursuing
the Ph.D. degree with the Institut National Polytech-
nique de Grenoble, Grenoble, France.
He joined Commisariat à l’Énergie Atomique–
Laboratoire d’Électronique et de Technologie de
l’Information, Grenoble, in 2012.
Gabriel Molas was born in Paris, France, in 1979.
He received the B.S. and M.S. degrees in physics
engineering with a specialization in microelectron-
ics, and the Ph.D. degree in micro and nanoelec-
tronics from the Polytechnics Institute of Grenoble,
Grenoble, France, in 2001 and 2004, respectively,
with a thesis in few electron memories.
He joined Commisariat à l’Énergie Atomique–
Laboratoire d’Électronique et de Technologie de
l’Information, Grenoble, in 2004, as a Research
Engineer.
Philippe Blaise received the master’s degree
in applied mathematics from the École
Nationale Supérieure d’Informatique et de
Mathématiques Appliquées de Grenoble,
Grenoble, France, and the Ph.D. degree
in physics from Université Joseph Fourier,
Saint-Martin-d’Hères, France, in 1998.
He is currently a Permanent Researcher with
Commisariat à l’Énergie Atomique–Laboratoire
d’Électronique et de Tec hnologie de l’Information,
Grenoble.
Mathieu Bernard , photograph and biography not available at the time of
publication.
Anne Roule received the Engineer degree in micro-
electronics with a specialization in physics of
materials from the Institut National des Sciences
Appliquées de Lyon, Lyon, France.
She joined Commisariat à l’Énergie Atomique–
Laboratoire d’Électronique et de Technologie de
l’Information, Grenoble, France, in 2003, as a Scien-
tist, where she is currently responsible for advanced
materials development.
Gilles Le Carval , photograph and biography not available at the time of
publication.
Vincent Delaye , photograph and biography not available at the time of
publication.
Alain Toffoli received the M.S. degree in engineer-
ing from the Conservatoire National des Arts et
Metiers, Paris, France.
He is currently the Head of the Advanced and
Statistical Electrical Characterization Test Team with
the Electrical Characterization and Test Labora-
tory, Commisariat à l’Énergie Atomique–Laboratoire
d’Électronique et de Tec hnologie de l’Information,
Grenoble, France. His current research interests
include test methodologies.
Gérard Ghibaudo received the Ph.D. and State
Ph.D. degrees from the Institut de Microélectron-
ique, Électromagnétisme et Photonique (IMEP)–
Laboratoire d’Hyperfréquences et de Caractérisa-
tion (LAHC) Laboratory, MINATEC Center, Greno-
ble, France, in 1981 and 1984, respectively.
He is currently the Director of Research with
the Centre National de la Recherche Scientiﬁque,
IMEP–LAHC Laboratory, MINATEC Center. His
current research interests include MOS device
physics, ﬂuctuations, and low-frequency noise and
dielectric reliability.
Fabien Clermidy received the master’s degree, the
Ph.D. degree in engineering science, and the Super-
visor degree from the Grenoble Institute of Tech-
nology, Grenoble, France, in 1994, 1999, and 2011,
respectively.
He is currently the Head with the Digital Design
Laboratory, Commisariat à l’Énergie Atomique–
Laboratoire d’Électronique et de Technologie de
l’Information, where he is involved in multicore
architectures and design with a focus on emerging
technologies.
Barbara De Salvo received the Ph.D. degree in
microelectronics from the Polytechnic Institute of
Grenoble, Grenoble, France.
She has been a Scientist with Commisariat
à l’Énergie Atomique–La boratoire d’Électronique
et de Technologie de l’Information, Grenoble,
since 1999.
Luca Perniola received the Laurea degree in
nuclear engineering from the Politecnico di Milano,
Milan, Italy, in 2002, and the Ph.D. degree from
the University of Pisa, Pisa, Italy, and the Insti-
tut National Polytechni que de Grenoble, Grenoble,
France, in 2005.
He was with Commisariat à l’Énergie Atomique–
Laboratoire d’Électronique et de Technologie de
l’Information (LETI), Grenoble, in 2005. He has
been appointed as the Head of the Advanced Mem-
ory Laboratory w ith LETI, s ince 2013.
Authorized licensed use limited to: BEIJING INSTITUTE OF TECHNOLOGY. Downloaded on January 05,2026 at 01:22:27 UTC from IEEE Xplore.  Restrictions apply. 