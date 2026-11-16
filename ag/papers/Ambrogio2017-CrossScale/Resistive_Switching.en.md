<!-- PDF page 1; local proof page count -->

J Electroceram
DOI 10.1007/s10832-017-0093-y
Modeling resistive switching materials and devices across
scales
Stefano Ambrogio1 · Blanka Magyari-K ¨ope2 · Nicolas Onofrio3 ·
Md Mahbubul Islam4 · Dan Duncan2 · Y oshio Nishi2 · Alejandro Strachan4
Received: 17 January 2017 / Accepted: 22 May 2017
© Springer Science+Business Media New York 2017
Abstract Resistance switching devices based on electro-
chemical processes have attractive significant attention in
the field of nanoelectronics due to the possibility of switch-
ing in nanosecond timescales, miniaturization to tens of
nanometer and multi-bit storage. Their deceptively simple
structures (metal-insulator-metal stack) hide a set of com-
plex, coupled, processes that govern their operation, from
electrochemical reactions at interfaces, diffusion and aggre-
gation of ionic species, to electron and hole trapping and
Joule heating. A combination of experiments and model-
ing efforts are contributing to a fundamental understanding
of these devices, and progress towards a predictive under-
standing of their operation is opening the possibility for
the rational optimization. In this paper we review recent
progress in modeling resistive switching devices at multi-
ple scales; we briefly describe simulation tools appropriate
at each scale and the key insight that has been derived from
them. Starting with ab initio electronic structure simula-
tions that provide an understanding of the mechanisms of
/envelopebackAlejandro Strachan
strachan@purdue.edu
1 Dipartimento di Elettronica, Informazione e Bioingegneria,
Italian Universities Nanoelectronics Team, Politecnico di
Milano, Milano 20133, Italy
2 Department of Electrical Engineering, Stanford University,
Stanford, CA, 94305, USA
3 Department of Applied Physics, The Hong Kong Polytechnic
University, Hong Kong, SAR, China
4 School of Materials Engineering and Birck Nanotechnology
Center Purdue University, West Lafayette, IN, 47906, USA
operation of valence change devices pointing to the impor-
tance of the aggregation of oxygen vacancies in resistance
switching and how dopants affect performance. At slightly
larger scales we describe reactive molecular dynamics sim-
ulations of the operation of electrochemical metallization
cells. Here the dynamical simulations provide an atomic
picture of the mechanisms behind the electrochemical for-
mation and stabilization of conductive metallic filaments
that provide a low-resistance path for electronic conduction.
Kinetic Monte Carlo simulations are one step higher in the
multiscale ladder and enable larger scale simulations and
longer times, enabling, for example, the study of variabil-
ity in switching speed and resistance. Finally, we discuss
physics-based simulations that accurately capture subtleties
of device behavior and that can be incorporated in circuit
simulations.
Keywords Resistive switching · Valence change
memory · Electrochemical metallization cells · Reactive
molecular dynamics · EChemDID
1 Introduction
Nanoscale devices that can undergo reversible resistance
switching driven by electrochemical processes hold great
promise as memory and logic elements in future nanoelec-
tronics devices. They have the potential to contribute to the
extension of Moore’s law beyond the rapidly approaching if
not ongoing end of scaling and even enable neuromorphic
computing [1–3]. Resistance- and threshold-switching cells
are fascinating devices that can reversibly switch between
low-resistance state (LRS) and a high-resistance state (HRS)
via the application of an external electrochemical potential.


---

<!-- PDF page 2; local proof page count -->

J Electroceram
Their simple structure, consisting of two metallic elec-
trodes separated by a solid dielectric or electrolyte, seems
at odds with the wide range of I-V characteristics that can
be achieved by the appropriate choice of materials [4]: from
linear to non-linear bipolar and nonpolar resistance switch-
ing to threshold switching [5] (an abrupt but reversible
change in resistance), see Fig. 1. These devices can exhibit
switching in nanosecond timescales [6] and scaling to
approximately ten nanometers [7], where the distinction
between material and device becomes meaningless. Ultra-
fast switching, scalability to the few nanometer regime
and compatibility with CMOS processing make resistance
switching materials especially attractive for memory appli-
cations. These resistive random access memories (RRAM)
store information using electrical resistance as the state vari-
able as opposed to the charge stored in a capacitor currently
used in SRAM, DRAM and Flash [8].
Excellent review papers and perspectives on the electrical
characteristics, physics, chemistry and application of these
devices have been published in recent years, see for example
Refs. [8, 9]; this paper presents a review the state of the art
in predictive modeling of electrochemical RRAM devices
across scales. We cover recent progress in theoretical and
simulation work, from electronic structure calculations that
revealed the role of point defects in resistance switching in
oxides, and reactive molecular dynamics capable of describ-
ing electrochemical reactions with atomistic resolution, to
continuum modeling of the electronic characteristics of
the entire devices. A brief description of the physics and
operation of these devices is provided for completeness.
Resistive switching devices The wide range of charac-
teristics exhibited by this class of devices is an exam-
ple of emergent phenomena originating from a set of
underlying materials processes: electrochemical reactions
at electrolyte/cathode and electrolyte/anode interfaces,
Resistance Switching Threshold Switching
Fig. 1 Atomistic structure of pristine ECM ( top left ) and one in
the low-resistance with a conductive filament. Bottom panels: various
possible I-V characteristics of electrochemical cells
electrochemical dissolution and deposition of metallic ions,
transport and aggregation of ions and defects, injection,
trapping and migration of electrons and holes, changes in
composition, phase and valence change in the electrolytes,
Joule heating and electromigration. How these underlying
processes conspire to result in the rich set of behaviors
observed remains not fully understood [9]. Due to the vari-
ety of dominant mechanisms these devices are loosely clas-
sified into electrochemical metallization cells [ 10]( E C M ) ,
valence change memory [11] (VCM), atomic switches and
threshold switching materials [4].
VCM cells involve oxides such as HfO
x ,T i Ox ,T a Ox
[12–14] separating electro-inactive metallic electrodes, their
simplicity and materials widely used in the semiconduc-
tor industry make them particularly attractive. A pristine
device is initially FORMED via the application of a voltage
larger than the operating voltage to induce a soft break-
down of the oxide and a transition to the LRS, is called
the forming step and is done normally at a voltage V
FORM
significantly higher than the device normal operating volt-
ages. Once a device has been formed, it can be reversibly
switched between the LRS state (ON) via the so-called SET
process and switched back to the HRS (OFF) by a RESET
process [13, 15–17]. The low resistance state has been asso-
ciated with the formation of a conductive filament of oxygen
vacancies that forms and dissolves driven by the injection
and trapping of electrons from the electrodes [12].
ECM cells operate in a similar fashion, but involve a
metallic filament. These devices consist of an electro-active
electrode (typically Cu or Ag) separated by an electro-
inactive one (often W, Pt or Au) by an oxide (SiO
2 or Al2O3)
or electrolyte (chalcogenides GeS/GeSe or Ag2S) [10]. The
application of an external voltage of appropriate polarity
results in the dissolution of the electro-active metal into the
dielectric, the transport of ions towards the inactive elec-
trode and their aggregation. This process eventually results
in the formation of a conductive metallic filament that short-
circuits the electrodes leading to a significant reduction in
resistance. A compliance current limits the current in the
ON state and its value governs the resistance of the ON state.
Increasing compliance current leads to lower resistance
(attributed to the formation of thicker filaments), enabling
multilevel bit storage and neuromorphic computing.
While in some cases debate persists about the switch-
ing mechanisms, interfacial reactions, and the nature of the
conducting paths, recent experiments and theoretical work
are helping build a quantitative picture of the operating pro-
cesses. Electrochemical tests like cyclic voltammetry have
provided information about the redox reactions as well as
the mobility and charge state of the ions [ 18, 19]. The
mechanisms associated with filament growth have also been
investigated intensively. Electrical measurements show that
the resistance of the ON state is independent of device

---

<!-- PDF page 3; local proof page count -->

J Electroceram
cross-sectional area, which is interpreted as indicating that
a single filament dominates transport but the decrease in
resistance with increasing compliance current could indicate
multiple filaments under certain conditions [ 6]. Given the
complexity of the processes involved in filament formation,
imaging has played an important role and high-resolution
transmission electron microscopy imaging and scanning
tunneling microscopy have provided key insight into the
formation of conductive filaments. [ 14, 20–22] Filaments
at various stages of formation have been imaged, includ-
ing under in situ conditions [23] providing important insight
into the formation of the conductive filaments. Despite such
progress, these experimental techniques are not without lim-
itations and the interpretation of the results often require
complex models; these limitations are particularly chal-
lenging when dealing with devices at their miniaturization
limit (few nanometers) and operating at ultra-fast speeds
(nanoseconds). Thus, physics based models at electronic,
atomic and device scales are a key complement to experi-
ments and are contributing to a complete understanding of
these fascinating devices.
Predictive, multiscale modeling of materials The descrip-
tion of materials from first principles without empirical
parameters became possible with the development of quan-
tum mechanics in the first decades of the 1900s. Paul Dirac
famously (and perhaps optimistically) stated in paper pub-
lished in 1929: “The fundamental laws necessary for the
mathematical treatment of a large part of physics and the
whole of chemistry are thus completely known, and the dif-
ficulty lies only in the fact that application of these laws
leads to equations that are too complex to be solved” [24].
Connecting the fundamental physics that describe the inter-
actions between electrons and ions to materials properties
and device characteristics has proved to be a gargantuan
endeavor and required significant developments in theory
and computation, from the development of efficient elec-
tronic structure methods, and potentials and algorithms for
large-scale molecular dynamics (MD) simulations to finite
elements. Bridging results across scales via multiscale mod-
eling is the avenue of choice today to realize Dirac’s vision,
see, for example. de Pablo and Curtin [25]. The field has
attracted significant attention in the last few decades as
compute power and advances in methods enabled quantita-
tive predictions for realistic materials. Multiscale modeling
often starts with a quantum mechanical description of elec-
tronic structure of the material of interest; methods like
density functional theory (DFT) are generally applicable
to any material of interest and provides results with the
accuracy required for many (but not all) applications [ 26].
Unfortunately, the computational intensity of these methods
preclude the simulation of large systems and long timescales
and a hierarchy of methods are often used to connect
ab initio predictions with materials and device performance.
The first step in this multiscale ladder is often MD simula-
tions with interatomic potentials developed from electronic
structure calculations; here electrons are averaged out and
the time evolution of the atoms in the material are solved
using classical equations of motion. Multi-million atom
MD simulations are common today and billion-simulations
are possible [27]. However, even using the largest super-
computers MD is unable to achieve the length scales and,
more often, the timescales required to describe many mate-
rials and device processes. Thus, information from ab initio
and MD simulations are often passed to continuum mod-
els capable of connecting fundamental processes to the
device world. Multiscale models have been developed to
describe, among other processes, mechanical [28–31], elec-
tronic [32, 33], chemical [34] properties of materials. In
this review paper we cover recent progress on modeling
resistive switching materials at the various scales mentioned
above and describe the insight, power and limitations of
each approach. Our aim is not to provide a comprehensive
review of all work in the field but to provide a perspective
of the current state of the art and challenges.
Manuscript organization We begin the paper with recent
advances in DFT studies of VCM, Section 2. These calcu-
lations have provided key insight into the nature, formation
and dissolution of the conductive filament and Section 3
discusses how doping can be used to engineer the response
of these devices. Recent advances in reactive MD simu-
lations enable the simulation of electrochemical processes
with time and spatial scales that match those of RRAM
devices at their scaling limit and operating at nanosecond
timescales. Section 4 describes reactive MD simulations and
Section 5 their
 application to ECM cells. Kinetic Monte
Carlo simulations enable achieving length and timescales
beyond those of MD while maintaining some description
of atomic processes and their variability, these methods are
discussed in Section 6. Since the early works, many models
have been developed to explain the resistive switching prop-
erties and the dependence from operational parameters such
as the external applied voltage or the maximum allowed cur-
rent in the device, also called as compliance current, or I
C .
Section 7 reviews a selection of the many physical models
of resistive switching, including both analytical and numer-
ical models and the observation of random telegraph noise
in Section 8.
2 Density functional theory: methods
An important class of RRAM devices, valence change
devices based on transition metal-oxide resistance-change
random access memory, and hafnium oxide (HfO 2)i n

---

<!-- PDF page 4; local proof page count -->

J Electroceram
particular, have shown great potential for applications in
next generation of nonvolatile memory and neuromor-
phic computing [35–40]. Along rapid developments on
the experimental side, quantum mechanical based simula-
tions [ 13, 15, 16, 41–51] have provided insight regarding
the nature of the conductive channel, the formation energy
of intrinsic defects and impurity ions, the drift and diffusion
mobilities of ions in a defective matrix, the charge trap-
ping effects near impurities and defects, and the electronic
structure of doped oxides.
In this section we briefly introduce density functional
theory, the theory underlying the calculation of electronic
structure of these materials, and how it is applied to char-
acterize the defects involved in VCM resistance switching.
Section 3 then describes the application of the technique to
understand the effect of doping on resistance switching.
2.1 Density functional theory: first principles
simulations of valence change memory devices
In density functional theory [52], the electronic structure is
often obtained by solving the Kohn-Sham equations [ 53].
They replace the full many body electron system with one
of particles that interact with each other via a mean field
approximation, and are subjected to an effective external
potential; this mapping is, in principle, exact for the ground
state density and energy. The Kohn-Sham equations are
expressed as:
(
− ℏ2
2me
∇2 + Veff(r)
)
/Psi1i(r) = ϵi/Psi1i(r), n(r) =
∑
i
|/Psi1i(r)|2, (1)
where Veff (r) is the effective potential that includes con-
tributions from the external potential due to the nuclei,
the electron-electron Coulomb potential and the exchange-
correlation potential, /Psi1i is the wavefunction of electron i
and ϵi is the energy of electron i.
The total energy of a system E[n] is expressed as a
functional of the charge density n(r)
E[n]= Ts[n]+
∫
drVext (r)n(r) +EHartree [n]+EII +EXC [n], (2)
where Ts is the Kohn-Sham kinetic energy of a system with
density n(r) assuming there are no electron-electron interac-
tions, Vext (r) is the external potential due to the nuclei and
any other external fields, EHartree is the Coulomb interac-
tion energy of the interacting electron density n(r) known
as the Hartree energy, EII is the interaction energy between
nuclei, and EXC is the exchange-correlation energy.
With this formalism, the ground state energy of a sys-
tem E[n] and its corresponding electron density n(r) can
be calculated if suitable approximations for the exchange-
correlation energy E
XC [n] and the external potential Vext (r)
can be found. Two early approximations for the exchange-
correlation energy EXC [n] were the local density approxi-
mation (LDA) and the generalized gradient approximation
(GGA). Developed later, the LDA + U approach is based
on LDA, but incorporates empirical on-site Coulomb cor-
rections. The simplified, rotationally invariant equation for
the LDA + U total energy is given by [54]:
E
LDA+U = ELDA + U − J
2
∑
σ
( ∑
m1
nσ
m1,m1
−
∑
m1,m2
ˆnσ
m1,m2ˆnσ
m
2,m1
)
(3)
Where E LDA and E LDA+U are the total energies cal-
culated using LDA and LDA + U, U and J are the
correction factors, m and σ are the magnetic and spin quan-
tum numbers, and n is the on-site occupancy matrix. Using
an optimal combination of orbital-specific correction factors
U and J , it is possible to perform calculations with a com-
putational demand typical for LDA or GGA, but reaching
desired accuracies on predicting band gaps and electronic
structures [43] comparable to hybrid functionals [55–58],
which include exact exchange.
2.2 Relative formation energy function for dopant
energetics in HfO
2
The defect types analyzed in this study are oxygen vacancies
(VO), cation substitutional dopants (D Hf), and interstitial
dopants (DI). The formation reaction in Kr ¨oger-Vink nota-
tion for a neutral cation substitutional dopant D on the Hf
site can be written as:
nHf
x
Hf + Dx
bulk
⇔ (n − 1)Hfx
Hf
+ Dx
Hf
+ Hfx
bulk
(4)
In order to facilitate comparison between different defect
sites, in the following the equations are rewritten in terms of
the chemical potential of oxygen, μO = 1
2O2(g):
For oxygen vacancy:
VO : Ox
O ⇔ Vx
O
+ 1
2O2(g) (5)
For substitutional dopant D on the Hf site:
DHf : (n − 1)Hfx
Hf + 2(n − 1)Ox
O
+ Dx
bulk
+O2(g) ⇔ (n − 1)Hfx
Hf
+ 2nOx
O
+ Dx
Hf
(6)

---

<!-- PDF page 5; local proof page count -->

J Electroceram
For interstitial dopant D:
DI : Hfx
Hf + 2Ox
O
+ Dx
bulk
⇔ Hfx
Hf
+ 2Ox
O
+ Dx
i
(7)
where n is the number of formula units in the supercell (e.g.
one unit cell of m-HfO2 contains four formula units). Thus,
the corresponding defect formation energy equations are:
VO : EForm (VO) = E(VO) −
(
E0 − 1
2 E(O2)
)
(8)
DHf : EForm (DHf) = E(DHf) −
((
n−1
n
)
E0 + E(D) + E(O2)
)
(9)
DI : EForm (DI) = E(DI) − (E0 + E(D)) (10)
Where EForm (x) is the formation energy of defect
x,E (VO)a n dE(Dx) are the total energies of the defective
HfO2 supercells containing defects V O and Dx, E0 is the
total energy of the pristine supercell, and E(O2)a n dE (D)
are the chemical potentials of O 2 and defect D. For the
chemical potential of O2 we used half of the binding energy
of an O2 molecule −5.375 eV . The advantage of using the
equations above for the EForm is that the relative formation
energy of a dopant at interstitial versus substitutional sites,
Er , is obtained effectively and the dopant chemical potential
dependence is eliminated.
Er(DI, DHf) = EForm (DI) − EForm (DHf)
= E(DI) − E(DHf) + E(O2) − 1
nE0 (11)
Er (DI,Dx) > 0 eV indicates dopant D favors the interstitial
site I over the substitutional site x,a n dEr (DI,Dx) > 0e V
indicates the substitutional site is favored.
2.3 Computational setup
Spin-polarized density functional theory calculations using
the Vienna ab initio Simulation Package (V ASP) were
carried out on monoclinic HfO 2 (m-HfO2, space group
P21/c) [59]. Projector augmented-wave pseudopoten-
tials [60], periodic boundary conditions, and an energy cut-
off of 367.38 eV were employed. K-point integration used
2x2x2 gamma-centered Monkhorst-Pack grids, and all ions
w e r er e l a x e dt oa ne n e r g yc o n v e r g e n c eo f1 0
−3 eV/atom
and forces less than 0.01 eV/ ˚A per ion. A unit cell of m-
HfO2 contains 12 atoms: 4 sevenfold-coordinated Hf, 4
threefold-coordinated (3C) O, and 4 fourfold-coordinated
(4C) O atoms. We used up to 3x3x3 supercells containing
324 atoms in our calculations. The optimum combination
of U parameters for monoclinic HfO
2 was found to be
Ud = 6.6e Va n dUp = 9.5e V[ 48, 49].
3 Role of doping on VCM switching
3.1 Conductive filaments and CVM general
characteristics
RRAM devices exhibit nonvolatile resistive switching under
an applied bias, they can reversibly switch between a high-
resistance state (HRS) and a low-resistance state (LRS),
representing the OFF and ON states. Switching in VCM
devices is hypothesized to be due to the formation and rup-
ture of a conductive filament through the oxide as depicted
in Fig. 2. An RRAM device that has never been resistively
switched is referred to as being in the pristine state. The first
switching step, which induces a soft breakdown of the oxide
and a transition to the LRS, is called the forming step and
is done normally at a voltage V
FORM significantly higher
than the device normal operating voltages. Once a device
has been formed, it can now be reversibly switched between
the ON and OFF states through the act of switching from
the HRS to the LRS referred to as the SET process, and
switched back from LRS to HRS by a RESET process [ 13,
15–17].
According to the prevailing theory of VCM operation
s h o w ni nF i g .2, the soft breakdown of the oxide is due to
the formation of a conductive filament of oxygen vacancies
through the oxide. The formation of the conductive filament
thus corresponds to a local reduction of the oxide (SET),
while the rupture of the filament corresponds to reoxida-
tion of the filamentary region (RESET) [ 13, 41]. The ratio
between the resistance in the HRS and that in the LRS
is referred to as the ON/OFF ratio, or alternatively as the
switching window. Generally, it is advantageous for a device
to have as large a switching window as possible. In contrast,
ECM devices operate via the formation and dissolution of
metallic filaments and ab initio simulations have been used
to address some of the of its characteristics features [61].
Our discussion here is focused on the bulk proper-
ties altered by the dopants, for a complete picture inter-
face effects should also be taken into account as recently
depicted by O’Hara [62]f o rH f / H f O
2 oxygen interface
diffusion. Moreover, the equilibrium charge state of the
dopants and vacancies will also affect the relative stability
as has been shown that oxygen vacancy charge states dra-
matically alter during the SET and RESET operation [ 16,
46, 50, 51].
3.2 Substitutional vs interstitial stability of dopants
in monoclinic HfO
2
For tuning device characteristics one possible option is to
introduce preferable dopants either during deposition or
by implantation. Doped RRAM devices have been shown

---

<!-- PDF page 6; local proof page count -->

J Electroceram
Fig. 2 Schematics of VCM
RRAM operation
to induce changes to the forming voltage, SET voltage,
power dissipation, ON/OFF ratio, variability, retention, and
endurance [ 44, 63–73]. Moreover, dopants reducing the
oxygen vacancy (VO) formation energy in their immediate
vicinity were linked to a strong impact on device behav-
ior [44, 64, 67, 73]. The preferred lattice site of dopant ions
will be investigated in detail to shed light on the implications
as substitutional dopants can contribute to V
O filament for-
mation and enhance VCM-type switching, while interstitial
dopants may form competing dopant conductive filaments,
which are not useful for VCM but are useful for ECM
cells (also called conductive bridge random access memory,
CBRAM) [74–76]. The relative energy of dopants in substi-
tional and interstitial ( Er ) based on Eq. 11 shows a strong
dependence on the dopant ion and valence. The more isova-
lent a dopant ion is with the species it is replacing, the more
stable the dopant will be on the substitutional site. And con-
versely, the more heterovalent ions tend to be more stable on
interstitial defect sites. From Fig.3, we observe a clear peri-
odicity in the relative formation energy of cation dopants,
pointing to isovalency with Hf as being a critical factor in
their relative stabilities. These calculations agree with exist-
ing experimental and other theoretical evidence [65, 67, 69,
70]: Al, Si, and Zr were reported to favor cation substitu-
tional sites in HfO
2, while interstitial sites are favored by
Ag, B, Cu, and Ni. On the other hand H, F, and N are all
known to favor anion substitutional sites.
3.2.1 V acancy formation energy in monoclinic HfO
2
with cation dopants
Oxygen vacancies in HfO 2 can be either in the 3C or 4C
coordinated configuration, with the 3C coordination being
slightly more energetically favorable. Figure 4 shows the
supercell of HfO 2 containing a single oxygen vacancy in
the 3C configuration. Next, to simulate the formation of
filaments in VCM devices a column of oxygen vacancies
was removed in a supercell HfO
2.I nF i g .5 we show how
this filament structure is modified with the introduction
of increasing concentration of titanium ions. To assess the
electronic charge transfer effects with the introduction of
dopants near the filament, the partial charge distribution
of the electronic states associated with the filament and
dopants are plotted in Fig. 6. This quantity is obtained by
calculating the spatial distribution of electronic states with
energies between the valence and conduction bands of m-
HfO
2, i.e. through defect states induced in the band gap. The
filament formation is clearly visible as a conductive chan-
nel formed by delocalized electrons from the metal ions in
the neighborhood of interacting oxygen vacancies. With the
introduction of Ti dopants a partial electronic localization is
observed in the immediate vicinity of Ti.
The ultimate goal of VCM doping is to tune the ener-
getics of oxygen vacancy formation within the device. In
this section, the dopant effect on V
O formation energy was
calculated for both 3C and 4C oxygen vacancies in HfO 2
according to Eq.8. As seen in Fig.7, the presence of dopants
generally reduce the formation energy of nearby V O ,h o w -
ever the magnitude of the change strongly depends upon
the dopant species. It was found that the effect of cation
dopants on V O formation energy depended heavily upon
their valence relative to Hf, with strongly n-type or p-type
dopants causing the greatest reductions in V O formation
energy. On the opposite, the effect of anions on VO forma-
tion is generally weaker and is less affected by the particular
dopant species.
Relative to their localized single vacancy effects, the
energetics and charge stability of oxygen vacancies when

---

<!-- PDF page 7; local proof page count -->

J Electroceram
Fig. 3 Relative energy for
substitutional and interstitial site
preference of dopants in
monoclinic HfO
2
they form extended defects like clusters or filaments, can
be significantly altered. This has been shown to take place
due to a number of factors including the Coulomb repulsion
between the vacancies, constraints on lattice relaxation, and
delocalization of unpaired Hf 5d electrons across neighbor-
ing Hf ions. Therefore, it is informative to inspect the impact
of dopants on the filament energetics. From Fig. 7 a similar
trend in the oxygen vacancy formation energy dependence
on dopants species is noted relative to their effect on sin-
gle vacancies. The overall effect is slightly weaker though,
because the formation energy is being spatially averaged
over the length of the entire filament, whereas the dopant
induces a strong local effect.
Fig. 4 Supercell structure of monoclinic HfO 2 with one oxygen
vacancy in the 3C position
Based on this data, implications on device performances
can be derived, which suggest that dopants presence in a
VCM switching oxide may induce preferential sites for the
formation of conductive filaments and therefore could be
beneficial to improve on the uniformity of devices. During
the forming step, the conductive channel nucleation can now
happen around dopants, since the formation energy of oxy-
gen vacancies near dopants is reduced compared to Hf-sites.
Moreover, if the same current compliance is used for un-
doped and doped devices then a more robust filament can
be achieved with dopants, which ultimately may improve
on the variability, retention and endurance characteristics of
VCM devices.
4 Reactive molecular dynamics simulations
of electrochemistry
MD simulations with reactive interatomic potentials enable
multi-million atoms simulations and the description of
Fig. 5 Oxygen vacancy filament in monoclinic HfO 2 with one and
two dopants

---

<!-- PDF page 8; local proof page count -->

J Electroceram
Fig. 6 Partial charge density of gap states for a monoclinic HfO2 with
an oxygen vacancy filament and two Ti dopant ions
complex chemistry, including combustion, decomposition
of explosives [77, 78], growth of carbon nanotubes [ 79],
among others [80]. The ability of reactive force fields like
REBO [81], COMB [82] and ReaxFF [83] to describe
the dissociation and formation of bonds lies in two key
concepts: i) the use of partial bond orders [ 84, 85]t o
describe covalent interactions, and ii) environment depen-
dent partial atomic charges that are computed based on
the atomic positions at each step of the simulation [ 86,
87]. This section introduces the electrochemical dynamics
with implicit degrees of freedom (EChemDID) method; an
extension of reactive MD simulations to describe electro-
chemical reactions driven by the application of an external
voltage between electrodes and, the generalization of the
reactive force field ReaxFF to describe discreet electron pro-
cesses. Section 5 describes the application of EChemDID to
simulate the dynamics of resistance switching in ECM cells.
4.1 Reactive MD simulations of electrochemistry:
EChemDID
The discussion of EChemDID here is based on the ReaxFF
force field, but the method is applicable to any reactive
potential (for example, a similar extension to COMB is
described in Ref. [ 88]). ReaxFF was chosen since it has
been successfully applied to model a wide range of materi-
als [80] including metals [89], semiconductors [90], organic
molecules [91] and interfaces [92].
Interatomic potentials describe the total energy of the
system as a function of atomic positions from which,
atomic forces are obtained to solve Newton’s equations of
motion in MD simulations. In ReaxFF, the total energy is
decomposed as the sum of many-body covalent interactions
(including bond stretch, bond angles, torsions) and pairwise
non-bonded interactions including van der Waals and elec-
trostatics based on environment-dependent partial atomic
charges. In the charge equilibration formalism (QEq) [ 86,
87], the total electrostatic energy is written as:
E ({R
i}, {qi}) =
∑
i
(
χ0
i qi + 1
2Hiq2
i
)
+
∑
i<j
qiqj J(Rij)
(12)
where χ0
i and Hi represent the electronegativity and hard-
ness of each atoms, respectively and J(Rij) is the shielded
(to avoid close range over-contributions) Coulomb function.
Interatomic charges {qi} are calculated self-consistently,
at every step of the simulation, by equating the elec-
tronic chemical potential {∂E/∂q
i} on each atom, assuming
charge neutrality.
Fig. 7 Formation energy of
single oxygen vacancies near
dopants and clustered oxygen
vacancies with one and two
dopants in monoclinic HfO
2


---

<!-- PDF page 9; local proof page count -->

J Electroceram
The effect of an electrochemical potential difference /Phi10
applied between two electrodes can be modelled by chang-
ing the electronegativity of the atoms on one electrode to
χ0
i → χ0
i + /Phi10/2a n dt oχ0
i → χ0
i − /Phi10/2 on the other. As
can be seen from Eq. 12, the external potential imposes an
energy difference of /Phi10e between electrons in either elec-
trode. The challenge in describing electrochemical reactions
lies in the dynamical evolution of the electrodes due to dis-
solution and deposition of atoms, consequence of Redox
reactions. Thus, the composition of each electrode evolves
in time; an atom in the electro-positive electrode of an ECM
cell may dissolve into the dielectric, migrate and merge to
the inactive electrode. In the original description [ 93], a
distance-based cluster analysis was performed to identify
the atoms belonging to each electrode and the electronega-
tivities were adjusted accordingly. This approach was later
generalized [ 94] by considering the local electrochemical
potential on atom i, /Phi1
i(t) as a dynamical variable whose
time evolution depends on the electrochemical potential of
nearby atoms. As a consequence, the atomic electroneg-
ativities used for charge equilibration also become time
dependent quantities χ
i(t)→ χ0
i + /Phi1i(t). Onofrio et al.
proposed a diffusive equation to describe the equilibra-
tion of electrochemical potential within connected metallic
regions:
˙/Phi1i(t) = k∇2/Phi1i(t) (13)
where k is an effective diffusivity for the electrochemical
potential and ∇2 is the Laplacian operator. In this approach,
the ultrafast equilibration of voltage, described rigorously
by Maxwell’s equations, is described via fictitious, but com-
putationally convenient, diffusion equation. The key to the
EChemDID approach is that the electrochemical potential
equilibration remains faster than the atomic processes. An
important by-product of using a diffusive equation is the
possibility to estimate electrical currents when a potential
difference is applied across a metallic system and one can
choose the effective diffusivity k to match the resistivity of
the electrode’s materials. A complete description of EChem-
DID is presented in Ref. [94] including a series of validation
tests and applications.
Practically, Eq. 13 is solved numerically at each step
of the simulation and the updated atomic electronegativity
χ
i(t) is used to perform the charge equilibration Eq. 12.
EChemDID is implemented as a LAMMPS user package
which can be downloaded from the authors’ webpage [95].
4.2 Explicit electrons in reactive MD simulations:
eReaxFF
The charge equilibration method discussed above describes
electron transfer between atoms, but is incapable of describing
discreet electron processes like the trapping of an electron
or a hole by a defect. As discussed above, these processes
are believed to critical in the operation of VCM and in
ECM cells where a dissolved Cu atom can be in the 1+ or
2+ states. An atomistic modeling for a proper theoretical
description of physical systems that involve explicit elec-
tron is also critical in rechargeable battery interfaces and
ferro-/piezoelectric materials. Ab initio electronic structure
methods are well established to provide the most accurate
description of such systems and processes. However, the
physical length and timescale required to model electron
flow driven processes are often prohibitive. The classical
reactive force field, such as ReaxFF, treatment of charge
transfer is inadequate to describe discreete electron transfer
associated reactive events [ 80, 96]. Recently, FF meth-
ods that include aspects of explicit valence electrons have
been developed, such as the electron force field (eFF) [ 97]
and the LEWIS [ 98] force-field. But, their capability to
describe complex reactions or materials have not been
demonstrated yet. The split-charge equilibration (SQE) has
also been extended to describe integer charge transfer in
reactive dynamics [ 99–101], but not applied to complex
chemical reactions yet. Islam et al. [ 102] recently extended
the ReaxFF method to incorporate an explicit electronic
DOF. In the eReaxFF development, basic framework and
functional forms of the ReaxFF were retained to ensure
transferability of the existing ReaxFF parameters to the
eReaxFF [102].
Theory A limited pseudoclassical electron/hole DOF has
been introduced in the eReaxFF method as a complement
to the implicit treatment of electrons in the ReaxFF method.
New energy functionals were designed to compute the elec-
trostatic interactions for explicit electrons and core charges.
The explicit electron- and hole- like particles carry neg-
ative (−e) and positive (+e) charges, respectively. The
pairwise electrostatic interaction between the electron and
core-charge was described as [103]
E
nucl(i)−elec(j) =− 1
4πεo
β
∑
i,j
Zi
Rij
erf(
√
2αRij) (14)
where Zi is the nuclear charge, Rij is the distance between
the electron and nucleus, α (Gaussian exponent), and β
are constants that depends on the atom type. The core-
charge is the charge corresponding to the atomic number
of an atom. In addition to the electron-electron interactions
as Coulomb point charges, short-range Gaussian repulsion
functions were also employed. In contrast to the ReaxFF
method, the explicit electron/hole results in variable valen-
cies and number of lone electron pairs of host atoms. The
variable atom valency was determined using an electron-
nuclear distance-dependent function; that also allows an

---

<!-- PDF page 10; local proof page count -->

J Electroceram
electron to virtually split itself among its neighboring atoms
as such resembles a partial delocalization. To accommodate
the electron-induced variable atom valency on the bonded
interactions, modified over/undercoordination energy func-
tional were used. The implementation of the charge-valency
coupling enabled eReaxFF method to describe chemical
bonding involving charged species more accurately via the
appropriate change in valence when calculating the degree
of over or under coordination of the host atoms. The ACKS2
charge calculation scheme [104] was implemented to cir-
cumvent unrealistic long range charge transfer issues of the
Electronegativity Equalization Method (EEM) [86].
The capability of the eReaxFF method has been demon-
strated through capturing electron affinities (EAs) of a
wide range of chemical species [ 102]. Figure 8 shows
that eReaxFF method qualitatively reproduces experimen-
tal [105, 106]a n dab initio values of EAs of all the saturated,
unsaturated, and radical species considered in the study,
which is a notable improvement in contrast to the corre-
sponding ReaxFF predictions. The authors employed the
developed force field to study electron transfer dynamics in
model hydrocarbon radicals, and good agreement with the
Ehrenfest dynamics simulations was established.
In Li-ion batteries (LIBs), during the first charge cycle,
reductive decomposition of the electrolyte molecules gener-
ates a thin passivation layer, which is known as solid elec-
trolyte interphase (SEI). Atomistic modeling of such reduc-
tion reactions requires treatment of explicit electrons. The
eReaxFF method enabled investigation of the SEI forma-
tion mechanisms via its explicit electronic DOF. In a recent
study, Islam et al. [107] applied the eReaxFF to describe
all the major reduction reaction pathways of SEI forma-
tion with ethylene carbonate (EC) molecule. The EC/Li
chemistry initiates through oxidation of lithium followed
by the reduction of EC. The EC
− /Li+ complex undergoes
EC ring opening reactions and produces o-EC− /Li+ radical.
The eReaxFF predicted energetics of the radical forma-
tion and termination reactions were in good agreement with
literature data. In dynamic simulations of EC/Li systems,
di-lithium butyl dicarbonate (Li
2BDC), and di-lithium ethyl
dicarbonate (Li2EDC) formation were observed, which are
considered as dominant components of anode-side SEI.
Figure 9 represents simulation snapshots of initial, gener-
ated o-EC
− /Li+ radicals, and Li 2BDC, Li2EDC formation
as observed during 600K NVT MD simulations. Overall, the
eReaxFF predicted EC/Li chemistries were in better agree-
ment with literature data compared to other ReaxFF studies
on LIBs electrolytes [96, 108]. The eReaxFF method with
its ability to adequately describe redox reaction is predicted
to be a powerful method for understanding the dynamics of
explicit electrons in physical and chemical systems.
5 Atomistic modeling of resistance switching
in ECM cells
Since ECM cells have been demonstrated to operate at ultra-
fast (nanosecond) timescale [6] and can be scaled down to
the nanoscale [7], EChemDID provides the unique tool
to study resistive switching with unprecedented spatial and
temporal resolutions. EChemDID was used to study Cu/a-
SiO
2 ECM cells; a common materials combination that
has received particular attention in the field because of its
compatibility with CMOS fabrication [ 109]. Details on the
atomistic models and annealing procedure are described
in Refs. [93, 110]. Here we describe the key atomistic
steps identified for stable resistance switching, discuss the
kinetics and stability of the nanoscale filament and show
how EChemDID can be extended and applied across the
scales of atomistic simulations.
Fig. 8 Electron affinities of
various species as calculated
using the eReaxFF method and
comparison with the
experimental [105, 106]a n d
ab initio data. Adapted from
[102]
-100
-80
-60
-40
-20
0
20
40
60
80
Electron Affinities (kcal/mol)

---

<!-- PDF page 11; local proof page count -->

J Electroceram
(a) (b) (c)
Fig. 9 Snapshot of the simulation cell at (a)t= 0p s ,(b) the generated
o-EC− /Li+ radicals are highlighted, t = 25 ps. EC and Li which are
not participated in the electron transfer event are displayed as line, (c)
The formation of Li2BDC, Li2EDC and C2H4 gas as observed in the
MD simulations. Color scheme: cyan: carbon, white: hydrogen, red:
oxygen, purple:L i+,l a r g eblue sphere: electron. Adapted from [107]
Figure 10 shows snapshots of the atomic structure of an
ECM cell operation corresponding to forming, reset and set
sequences. The switching mechanism can be summarized as
follow. When a voltage of the appropriate polarity is applied
between metallic electrodes, Cu atoms dissolve inside the
solid electrolyte, often aggregate into small clusters and are
electric-field driven toward the inactive electrode, where
they undergo successive reduction reactions. This process
leads to a nanoscale filament growing from the inactive elec-
trode toward the active one until a stable metallic bridge is
formed. When the voltage is reversed, the filament breaks
at the inactive electrode leaving a partial filament connected
to the active electrode. During a subsequent set voltage,
the partial filament at the active electrode rapidly dissolves
and another metallic bridge is established within shorter
timescale than that of the forming sequence. Details on the
mechanism are discussed in Ref. [93].
While EChemDID simulations, at sub-10 nm and ultra-
fast scales, indicate the formation of small Cu clusters and
their stabilization by reduction to be the rate-limiting step,
depending on materials and conditions, experimental stud-
ies show the existence of a critical nucleus size (ranging
from few to tens of atoms) corresponding to a nucleation
and growth process [22, 111, 112]. DFT calculations show
that aggregation of Cu atoms into clusters to be energetically
favorable in silica for all sizes without a critical nucleus
size. Thus, given the relatively low activation energy for dif-
fusion, clusterization is the rate limiting in such nanoscale
device [113]. We note that these DFT simulations consider
only neutral states under no bias and it is often observed
FORMING
RESET SET
FORMING RESET SET
0                         1                       2                         3                         4                           5                        6          (ns)
HR
LR
+8V -8V +8V
Fig. 10 Snapshots of the atomic structure during forming, reset and
set phases of a Cu/SiO2 ECM cell. The top scale represents the switch-
ing state of the cell as a function of the simulation time ranging from
high resistivity (HR) to low resistivity (LR). The colors represent
atomic charges ranging from −0.3e (red )t o +0.3e (blue). The solid
electrolyte has been hidden for clarity. Adapted from Ref. [93]

---

<!-- PDF page 12; local proof page count -->

J Electroceram
that nucleation is mainly rate limiting at low applied volt-
ages (up to 0.5 V) [18, 111, 112] and the voltage window
depends on the materials and operating conditions.
Simulations over an ensemble of statistically indepen-
dent ECM cells show that approximately 50% of the devices
can be formed in less than 5 ns under 8 V and that the dis-
crepancy in switching times is directly due to the atomic
variability of the amorphous structure of the electrolyte.
Interestingly, the dissolution and aggregation of Cu ions
inside silica does not generate point defects, ruling out
the formation of permanent nanoscale voids and channels
during the forming phase. The comparison of switching
timescales between different electrode geometries suggest
that nanoscale roughness leads to larger dissolution rate and
therefore shorter switching times [93, 110]. Finally, various
shapes of the nanoscale filament have been observed and
many simulations suggest that single atom chain bridges are
metastables at room temperature.
Although EChemDID does not describe electrons explic-
itly, it is possible to extend its application across the scales
by coupling the formalism with first principles and semi-
empirical methods. The evolution of the electronic density
of states (DoS), calculated from DFT, at various timestep
during switching of a small ECM cell with EChemDID
shows an increase of the number of Cu states at the
Fermi energy. The change in the electronic structure of
the electrolyte is directly responsible for transport in the
device [110, 114]. However, the high computational cost of
DFT restricts this approach to small ECM cells composed
of few hundreds of atoms. In another approach [ 115], the
(QEq) charges generated during an EChemDID simulation
were used to calculate the electrostatic potential described
by the Poisson equation and, an empirical tight-binding
(TB) Hamiltonian was reconstructed. Following this proce-
dure, the electronic structure of large ECM cells has been
resolved and accurate electronic currents have been com-
puted from the transmission spectrum. The TB study of
Cu/a-SiO
2 ECM cells shows a dramatic increase of the
transmission at the Fermi energy as soon as Cu atoms start
to dissolve and clusterize inside the electrolyte. A current
ratio of approximately 200 was calculated between pristine
and formed cells.
To conclude, EChemDID has enabled the first atom-
istic simulations of ECM cell operation and, the analysis
of switching in Cu/a-SiO
2 cells has provided a clear pic-
ture of the corresponding atomic mechanism. EChemDID
is not without limitation and, the description of atomistic
electrochemistry from a classical theory involves obviously
some approximations. An improved version of the reactive
force field including integer charge transfer and the explicit
treatment of electron which has been discussed in the pre-
vious section appears as one of the most appealing route to
improve EChemDID.
6 Kinetic Monte Carlo models
Timescales beyond tens of nanoseconds are beyond the
reach of MD simulations and exploring such regimes (for
example multiple switching events or slower switching)
requires coarser-grained approaches. One such method is
kinetic Monte Carlo (KMC) that describes the stochastic
evolution of a system given possible reactions and the asso-
ciated transition rates, which can be computed from ab
initio simulations. KMC has been applied to VCM cells and
describe the resistance switching through a random network
of resistances [116]o rt r a p s[117–119], each contributing to
conduction by trap-assisted tunneling. Also ECM systems
have been studied through KMC simulations [120–123].
Guan et al. [ 118] proposed a model to describe switch-
ing dynamics and variability. In this model, the oxide layer
between the electrodes is considered as oxygen-vacancy
(V
O )r i c h .VO behave as traps which assist electronic con-
duction through trap-assisted-tunneling, or TAT. Depending
on the concentration and distribution of V O , different con-
figurations are obtained, which lead to different high- and
low-resistance states, therefore giving an interpretation of
RRAM variability of the resistive states.
In this framework, the application of a voltage causes the
migration of oxygen ions, leading to generation and recom-
bination of V
O and thus causing current and resistance vari-
ability. The study is performed on a 2D geometry, as shown
in Fig. 11. Recently, 3D KMC models have also been pro-
posed [124, 125]. Figure 11(a) shows the area between the
biased electrode, on the right, and the grounded electrode,
on the left. The random points between them represent the
initial V
O before forming. Figure11(b) shows the calculated
forming current, while Fig. 11(c) shows the new V O gen-
erated by local fields, which arise when biasing the device.
The V
O generated create conductive filaments and perco-
lation paths where electronic conduction is favored. The
model is also capable of explaining the increasing resistance
states obtained by applying decreasing negative stop volt-
ages, namely, the maximum negative voltage applied to the
device. Figure 11(d) shows the calculated reset I-V curves,
which display a gradual reset, and corresponding increasing
high-resistance states for decreasing stop voltage from−2.5
Vt o−3 V . This dependence is interpreted by the formation
of a gap with a negligible concentration of V
O due to the
VO -oxygen recombination, since negative oxygen ions are
rejected by the negative electrode on the right. The formed
gap obstacles electrons motion, since there are no traps
aiding electron transport, thus rising the device resistance.
Current fluctuations during reset transition are also
explained by the model: in the gap, the electric field is
high, hence sometimes there is some V
O generation which
allows electrons to reach the left electrode, increasing the
current. However, since the recombination probability is

---

<!-- PDF page 13; local proof page count -->

J Electroceram
(a) (b) (c)
(d) (e) (f)
(g) (h) (i)
Fig. 11 RRAM switching model presented by Guan et al. [118].
(a) shows a 2D plane with the pre-forming oxide, with few VO .T h e
forming I-V curve (b) leads to the creation of percolation paths formed
by VO (c). Decreasing the stop voltage leads to increasing high resis-
tance levels (d) due to the increased tunneling gap. The tunnel length
increases for decreasing stop voltage (e–f ), reducing the probability
of electron transport through the gap. The model also accounts for set
transition (g). Higher I
C leads to the creation of more percolation paths
(h–i). Adapted from [118]
high because the oxygen concentration is elevated at the
left electrode, the generated V O will rapidly recombine.
In addition, the model can also account for set transition:
Fig. 11(g) shows two set I-V curves for different compliance
currents IC. Higher IC leads to the creation of more percola-
tion paths, as shown in Fig.11(h,i), thus lowering the device
resistance.
7 Continuum models of filamentary resistance
switching
7.1 Analytical models
Analytical models have a particular appeal due to the abil-
ity to describe complex physical processes with just a few
equations which can be easily solved. The main idea is to
describe only the relevant dependence of resistive switching
from the evolution of few parameters, such as CF diam-
eter [126–128], CF oxide gap [128–134]o rC Ft o pa n d
bottom reservoir [ 135]. Note that the three are geometri-
cal parameters, which allow for an immediate translation in
electrical resistance. The fast evaluation of resistance and
current dependence on voltage allows to implement such
models in circuit simulators, e.g., SPICE, thus enabling the
study of Monte-Carlo variability issues [136–138], memory
[1, 139, 140] logic [ 9] and novel neuromorphic networks
[141–143] with a high degree of accuracy, at a manageable
computational load.
I nt h ew o r kb yA m b r o g i oe ta l .[129], resistive switch-
ing is geometrically described as the growth of a CF with
an increasing diameter φ during set transition, namely the
switching from the high- to the low-resistance state. In the
same model, reset transition, namely the reverse switching

---

<!-- PDF page 14; local proof page count -->

J Electroceram
from low- to high-resistance state, is instead obtained by
means of an oxide gap growth of length /Delta1inside the CF,
which causes the rupture of the conductive path, leading
to an increased resistance. A schematical representation of
the reset transition is shown in Fig. 12. Figure 12(a) shows
a full set state, represented by a cylinder, which also cor-
responds to the initial condition. The CF is composed by
defects, for example oxygen vacancies V O , metal ions (Hf
in HfOx RRAM) or other impurities. The top and bottom
electrodes are fabricated in TiN, and the top electrode shows
a Ti layer, also called vacancy-rich layer, or, equivalently,
oxygen depletion layer, which acts as a defects reservoir.
The application of a negative voltage −V
A at the top elec-
trode leads to a migration of defects driven by the applied
electric field and accelerated by temperature. The migration
starts from the CF point at highest temperature, which cor-
responds to the center of the CF, Fig. 12(b), since top and
bottom electrodes act as heat sinks. This leads to a subse-
quent opening of an oxide gap of length/Delta1,F i g .12(c), which
increases the CF resistance. The growth rate of the gap can
be modeled by an Arrhenius process:
d/Delta1
dt = Ae
− EA
kTz1 , (15)
where A represents a pre-exponential factor and k is the
Boltzmann constant. The energy barrier E
A, which defects
need to overcome to migrate, is lowered by the applied
electric field:
E
A = EA0 − αqV, (16)
(a) (b) (c)
(d) (e) (f)
Fig. 12 Geometrical representation of the reset transition ( a-b-c)
which is caused by the gradual opening of an oxide gap, and of the set
transition (d-e-f), characterized by the growth of a sub-CF inside the
gap region. Adapted from [129]
and EA0 is the zero-field barrier, showing values around
1e V[ 144]. α is the barrier lowering factor, q is the electron
charge and V is the applied voltage on the gap region, which
drives the migration. Temperature T in Eq. 15 is calculated
by solving the steady-state Fourier equation:
kth
d2T
dz2 + J 2ρ = 0, (17)
which is consistently calculated in the three domains (top
CF stub, gap, bottom CF stub) by using the corresponding
CF and gap thermal conductivities k th [145] and resistivi-
ties ρ. J is the current density flowing into the CF region.
In Eq. 15, temperature is calculated at the bottom CF stub,
z1, since migration of defects starts here, hence z 1 repre-
sents the defect injection point. The final reset state shown
in Fig. 12(c) also acts as the initial state for the set transi-
tion, Fig. 12(d). During set transition, a sub-CF of diameter
φ nucleates and grows inside the gap region, Fig. 12(e),
eventually leading to a complete set state, Fig. 12(f). Since
defects are still field driven and temperature activated, the
corresponding φ growth rate can be still modeled as an
Arrhenius process:
dφ
dt = Ae
− EA
kTz2 . (18)
Parameters are identical to the reset transition, but tempera-
ture is now calculated at the z 2 edge, since now defects are
injected from the top CF stub towards the bottom electrode.
The model I-V curve result is shown in Fig. 13, obtaining
good quantitative results, particularly regarding switching
voltages V
set and V reset , corner voltage V C and reset
Fig. 13 Experimental and calculated with the analytical model in
[129] I-V curves. The model can capture the abrupt set transition and
the gradual reset transition. Adapted from [129]

---

<!-- PDF page 15; local proof page count -->

J Electroceram
current Ireset . Abrupt set transition is explained through a
positive thermal feedback: as the sub-CF grows, the cur-
rent flowing increases, leading to a temperature increase
which accelerates defects migration and the consequent CF
growth, thus leading to more current. Also gradual reset is
well captured. It is instead modeled through a negative ther-
mal feedback: the growth of an oxide barrier increases the
resistivity and cools the CF. Therefore, to activate a further
defect migration and continue the reset transition, voltage
must be increased to increase T again.
Another example of analytical model is the dynamic
‘hour glass’ model for resistive switching proposed by
Degraeve et al. [135]. In this model, set and reset switching
kinetics are caused by the migration of oxygen vacancies
V
O . The model schematizes the CF geometry as a V O top
reservoir (TR), a bottom reservoir (BR) and a thin bottle-
neck called constriction where conduction happens through
a quantum point contact, as shown in Fig. 14.
7.2 Numerical models
Numerical models represent an extension of the one-
dimensional analytical models previously shown, since they
solve a set of equations on a 2D or 3D domain [146–148].
From one side, the computational load increases due to a
more complex geometry and a denser space mesh, which
poses limitations in the use of numerical models for Monte
Carlo approaches, e.g., for variability and statistics. From
the other side, however, numerical models allow for a much
deeper understanding and description of the physical pro-
cess. The partial differential equations describing the system
physics, e.g., charge and defect transfer and heat transfer,
are consistently solved through a finite element method in a
2D or 3D geometry. Many numerical models have been pro-
posed, which include boundary effects, such as the impact
of electrodes [ 149] or space charges [ 150]. An example of
numerical model is the one which has been proposed by
Larentis et al. [146]. The switching mechanism is explained
as the motion of ionized defects, following drift-diffusion
equations. In this model, the electrostatic field is calculated
by solving the Poisson equation:
∇σ∇/Psi1= 0, (19)
where σ =1 /ρ is the electrical conductivity and /Psi1is the
local electrostatic potential. Note that no free space charge
has been considered, therefore the second term is zero. /Psi1is
linked to the electric field F through:
−∇ /Psi1= F. (20)
Heat transfer is taken into account by solving the steady-
state Fourier equation:
−∇ k
th∇T =| σ∇/Psi1|2. (21)
The right-hand side of Eq. 21 represents the local dissipated
power density, which is compared to the space variation
of the heat-flow caused by thermal conduction in the left-
hand side of Eq. 21. No thermal transient effects have been
considered since, for the CF region, the thermal constant
is much faster than any resistance switching effect [ 126].
The calculated temperature allows to obtain both diffusiv-
ity D and mobility μ of defects, which are related through
the Einstein relation and are thermally activated. All these
parameters are linked in the drift-diffusion equation which
describes the motion of defects during resistive switching:
∂n
D
∂t =∇ (D∇nD − μFnD), (22)
where nD represents the concentration of defects. Finally
both local kth and σ depend on nD.S i n c enD defects behave
like a dopant in an oxide layer, a high concentration of n D,
as in the CF, leads to high values of k th and σ. Instead,
depleted regions with a low nD concentration behave like an
oxide layer, with low kth and σ.
Figure 15 plots both experimental and calculated reset
I-V curves, showing a good agreement. Note that reset tran-
sition is shown for positive voltages, however this is due to a
different voltage sign convention, no physical changes arise
with respect to Fig. 13. The two experimental I-V curves
were obtained for different I
C, causing two initial low resis-
tance states of 0.4 and 1 k/Omega1. Increasing the voltage from 0
up to 0.4 V leads in a linear behavior since the cell shows
no migration of defects. At V = V
reset = 0.4 V , resistance
starts to increase due to the onset of migration activated by
temperature and field driven. Then, for increasing voltage,
resistance gradually increases up to a saturation of the resis-
tance at around 1 V . The model can closely capture these
Fig. 14 Initial device structure
(a), after forming (b)a n d
schematized geometry adopted
in the “hour glass”i model (c).
Adapted from [135]
(a) (b) (c)

---

<!-- PDF page 16; local proof page count -->

J Electroceram
Fig. 15 Measured and calculated reset I-V curves for different I C .
Figure 16 refers to points A, B, C and D. Adapted from [146]
dynamics also for different initial resistance states. Four dif-
ferent points are highlighted in the I-V curve of Fig. 15
corresponding to an initial resistance of 0.4 k /Omega1, namely
points A, B, C and D along the calculated curve. These
points are selected to study the evolution of reset transition.
Figure 16 shows the corresponding 2D maps calculated with
the numerical model related to the four points for (a) defect
concentration n
D in the CF region, (b) temperature T and
(c) potential /Psi1profile. At point A, the applied voltage is
slightly below reset voltage, therefore immediately before
the reset onset. The CF n
D concentration is almost contin-
uous and the temperature profile shows a parabolic shape
with a temperature maximum at the centre of the CF. The
corresponding potential /Psi1map reveals a uniform voltage
drop on the CF region. Immediately after reset, which is rep-
resented by point B, the CF is broken and a gap appears in
the n
D map, with a migration of defects from top to bot-
tom, as shown by the blue gap and red high density spot.
Note that the defect migration is opposite to Fig. 12 since
the geometry is rotated. The temperature shows an increase
in correspondence of the gap region, causing a thermal acti-
vation of migration. Also the potential starts to be localized
on the same region, which is consistent with a field-driven
defect migration. This localization is more pronounced for
higher voltages, which can be evidenced by maps at points
Ca n dD .
8 Random telegraph noise
One of the major concerns regarding RRAM reliability is
the ability to properly work at low current, e.g below 10
μA. In this regime, variability in switching time and resis-
tance increase [136] resulting in Random Telegraph Noise
Fig. 16 Calculated 2D maps of
defect concentration nD (a),
temperature T (b) and potential
/Psi1(c) for increasing reset states
A, B, C and D in Fig. 15.
Adapted from [146]
(a)
(b)
(c)

---

<!-- PDF page 17; local proof page count -->

J Electroceram
(a) (b)
Fig. 17 Experimental measurement of RTN for different I C ,( a).
When IC decreases, the relative amplitude of resistance fluctuation
/Delta1R=Roff -Ron increases. The measured and calculated relative ratio
/Delta1R/R as a function of R is shown in (b) for different materials,
revealing an initial slope for low resistance, followed by a saturation
at R > 100 k/Omega1. The inset shows a schematic representation of the CF
with a depleted region caused by the bistable defect. Adapted from
[151]
(RTN) and complicating the read operation. As an exam-
ple, Fig. 17(a) shows three different read resistances for
different compliance currents I C of 45 μA, 12 μAa n d
4 μA. As discussed above, decreasing I C results in an
increase in overall resistance but also in /Delta1R the differ-
ence between high R off and low R on resistance levels.
Figure 17(b) shows the corresponding normalized /Delta1Ra sa
function of the mean R of the device. Data from different
materials, HfO
2, NiO, and Cu-based CBRAM, show a lin-
ear increase for low resistance, followed by a saturation at
high resistance, above 100 k /Omega1. To tackle this dependence,
several models have been proposed [151–155]. The numeri-
cal model proposed by Ambrogio et al. [151], considers the
CF as a semiconductor-like material. Through a 3D numer-
ical simulator, electron transport equations are consistently
(a) (b)
Fig. 18 Calculated electron density nD maps inside the CF region, for
low-( a)a n dhigh-( b) resistance states. Adapted from [151]
solved to obtain the flowing current. In this model, the pres-
ence of a bistable trap near or inside the CF modifies the
device potential profile. In particular, the trap, which can be
a point defect, a dangling bond, or some impurity, can fluc-
tuate between a neutral and a negatively charged state. In its
neutral state, the defect does not affect electronic conduc-
tion; however, when it is negatively charge a section of the
CF is depleted, reducing the effective CF diameter available
for conduction and rising the device resistance. The calcu-
lated electron density n
D concentration maps are shown in
Fig. 18 for two different CF sizes, hence corresponding to
different IC. Figure 18(a) represents a relatively high I C :
the depleted region, represented in blue, partially depletes
the CF cross-section, giving rise to the linear slope in the
calculated results in Fig. 17(b). Figure 18(b) shows instead
al o w - I
C CF. In this latter case, the CF cross-section is
completely depleted, hence the difference between Ron and
Roff is high, giving a saturated regime for high resistance
in Fig. 17(b). In addition, Fig. 17(b) reports different calcu-
lated results for variable CF doping density N D, from 1018
to 1021 cm−3 . Negligible differences arise, meaning that the
size-dependence of RTN is mainly a geometrical effect.
9 Conclusions and outlook
Modeling and simulations at scales ranging from electrons
to those of the device have contributed significantly to
our understanding of the operation of resistance switch-
ing devices. From the characterization of the electronic and
atomic unit processes that result in resistance switching to
the interplay between heat generation and dissipation and
the growth and dissolution of a conductive filament at the
device level, ab initio simulations, molecular dynamics and
continuum modeling are today an indispensable partner to

---

<!-- PDF page 18; local proof page count -->

J Electroceram
experiments. Both simulations and experimental techniques
have limitations, involve approximations and assumptions,
and only their synergistic combination will result in a com-
plete understand of these devices and enable their rational
optimization.
Ab initio simulations, such as those described in Sections2
and 3 provide the most fundamental and broadly appli-
cable description of RRAM. Their explicit description of
electronic degrees of freedom enabled, for example, the
description of how charge trapping results in the formation
of conductive filaments in VCM. In addition, being applica-
ble to a wide range of materialsab initio methods have been
used to explore how doping affects the performance of these
devices. While very accurate, ab initio simulations can-
not capture the spatial and temporal scales involved in the
switching of realistic devices. Fortunately, MD simulations
with reactive interatomic potentials provide a less compu-
tational intense alternative and dynamical switching with
spatial and temporal scales matching those of devices at
the limit of scalability were recently demonstrated in ECM
cells, Section 5. These simulations simplify the description
of electronic degrees of freedom and often neglect discreet
electron events however, recent developments (Section 4)
show promise in the ability to capture such processes. One
step higher in the multiscale hierarchy, kinetic Monte Carlo
models, Section 6, have been developed to describe the
migration and aggregation of charge species. These models
can be parameterized from ab initio and MD simulations
and enable studying multiple switching and characterize the
intrinsic variability in the devices. Finally, numerical and
analytical continuum models, Section7 follow the evolution
of the conductive filament and temperature fields in the device
and can be used within circuit simulations in design experiments.
We are unaware of the systematic integration of the vari-
ous models and data across scales within a predictive mod-
eling framework for the prediction, design and optimization
of RRAM. We believe that the individual models are mature
enough for this integration and that such a modeling frame-
work together with appropriate experimental fabrication and
characterization can accelerate the development of RRAM
devices with improved properties.
Acknowledgments AS acknowledges support from the FAME Cen-
ter, one of six centers of STARnet, a Semiconductor Research Cor-
poration program sponsored by MARCO and DARPA. AS, NO and
MI acknowledge support for computational facilities from Purdue
University, the Center for Predictive Materials and Devices and the
Network for Computational Nanotechnology. BMK, DD and YN
acknowledges funding from the Stanford Non-V olatile Memory Tech-
nology Research Initiative Center and the calculations were carried out
using computational resources awarded on the NSF funded XSEDE
computational project and on the Carbon cluster at the Center for
Nanoscale Materials supported by the U. S. Department of Energy,
Office of Science, Office of Basic Energy Sciences, under Contract
No. DE-AC02-06CH11357.
References
1. R. Waser, M. Aono, Nanoionics-based resistive switching mem-
ories. Nat. Mater. 6(11), 833 (2007)
2. S.H. Jo, T. Chang, I. Ebong, B.B. Bhadviya, P. Mazumder, W.
Lu, Nanoscale memristor device as synapse in neuromorphic
systems. Nano Lett. 10(4), 1297 (2010)
3. M. Prezioso, F. Merrikh-Bayat, B.D. Hoskins, G.C. Adam, K.K.
Likharev, D.B. Strukov, Training and operation of an integrated
neuromorphic network based on metal-oxide memristors. Nature
521(7550), 61 (2015)
4. R. Waser, R. Dittmann, G. Staikov, K. Szot, Redox-based resis-
tive switching memories–nanoionic mechanisms, prospects, and
challenges. Adv. Mater. 21(25-26), 2632 (2009)
5. A. Padilla, G.W. Burr, R.S. Shenoy, K.V . Raman, D.S. Bethune,
R.M. Shelby, C.T. Rettner, J. Mohammad, K. Virwani, P.
Narayanan, et al., On the origin of steep–nonlinearity in mixed-
ionic-electronic-conduction-based access devices. IEEE Trans.
Electron Devices 62(3), 963 (2015)
6. Y .C. Yang, F. Pan, Q. Liu, M. Liu, F. Zeng, Fully room-
temperature-fabricated nonvolatile resistive memory for ultrafast
and high-density memory application. Nano Lett. 9(4), 1636
(2009)
7. V . Zhirnov, R. Meade, R.K. Cavin, G. Sandhu, Scaling limits of
resistive memories. Nanotechnology 22(25), 254027 (2011)
8. H.S.P. Wong, S. Salahuddin, Memory leads the way to better
computing. Nat. Nanotechnol. 10(3), 191 (2015)
9. J.J. Yang, D.B. Strukov, D.R. Stewart, Memristive devices for
computing. Nat. Nanotechnol. 8(1), 13 (2013)
10. I. Valov, R. Waser, J.R. Jameson, M.N. Kozicki, Electrochemical
metallization memories – fundamentals, applications, prospects.
Nanotechnology 22(25), 254003 (2011)
11. K. Szot, W. Speier, G. Bihlmayer, R. Waser, Switching the elec-
trical resistance of individual dislocations in single-crystalline
SrTiO
3. Nat. Mater. 5(4), 312 (2006)
12. J.J. Yang, F. Miao, M.D. Pickett, D.A.A. Ohlberg, D.R. Stew-
art, C.N. Lau, R.S. Williams, The mechanism of electroforming
of metal oxide memristive switches. Nanotechnology 20(21),
215201 (2009)
13. B. Magyari-K ¨ope, M. Tendulkar, S.G. Park, H.D. Lee, Y . Nishi,
Resistive switching mechanisms in random access memory
devices incorporating transition metal oxides: TiO 2,N i Oa n d
Pr0.7Ca0.3MnO3. Nanotechnology 22(25), 254029 (2011)
14. A. Wedig, M. Luebben, D.Y . Cho, M. Moors, K. Skaja, V .
Rana, T. Hasegawa, K.K. Adepalli, B. Yildiz, R. Waser, et al.,
Nanoscale cation motion in taox, HfOx and TiOx memristive
systems. Nat. Nanotechnol. 11(1), 67 (2016)
15. S.G. Park, B. Magyari-K ¨ope, Y . Nishi, Impact of oxygen vacancy
ordering on the formation of a conductive filament in for resis-
tive switching memory. IEEE Electron Device Lett. 32(2), 197
(2011)
16. K. Kamiya, M.Y . Yang, B. Magyari-K¨ope, M. Niwa, Y . Nishi, K.
Shiraishi, Physics in designing desirable reRAM stack structure
– atomistic recipes based on oxygen chemical potential con-
trol and charge injection/removal, inIEEE International Electron
Devices Meeting (IEDM) Technical Digest (IEEE, 2012), p. 2012
17. K.H. Xue, B. Traore, P. Blaise, L.R.C. Fonseca, E. Vianello, G.
Molas, B. De Salvo, G. Ghibaudo, B. Magyari-K ¨ope, Y . Nishi,
A combined ab initio and experimental study on the nature of
conductive filaments in resistive random access memory. IEEE
Trans. Electron Devices 61(5), 1394 (2014)
18. C. Schindler, G. Staikov, R. Waser, Electrode kinetics of Cu-
SiO
2-based resistive switching cells Overcoming the voltage-
time dilemma of electrochemical metallization memories. Appl.
Phys. Lett. 94(7), 2109 (2009)

---

<!-- PDF page 19; local proof page count -->

J Electroceram
19. S. Tappertzhofen, H. M ¨undelein, I. Valov, R. Waser, Nanoionic
transport and electrochemical reactions in resistively switching
silicon dioxide. Nanoscale 4(10), 3040 (2012)
20. C. Schindler, S.C.P. Thermadam, R. Waser, M.N. Kozicki, Bipo-
lar and unipolar resistive switching in Cu-doped SiO 2. IEEE
Trans. Electron Devices 54(10), 2762 (2007)
21. Y . Yang, P. Gao, S. Gaba, T. Chang, X. Pan, W. Lu, Observation
of conducting filament growth in nanoscale resistive memories.
Nat. Commun. 3, 732 (2012)
22. I. Valov, I. Sapezanskaia, A. Nayak, T. Tsuruoka, T. Bredow,
T. Hasegawa, G. Staikov, M. Aono, R. Waser, Atomically con-
trolled electrochemical nucleation at superionic solid electrolyte
surfaces. Nat. Mater. 11(6), 530 (2012)
23. W.A. Hubbard, A. Kerelsky, G. Jasmin, E.R. White, J. Lodico, M.
Mecklenburg, B.C. Regan, Nanofilament formation and regener-
ation during Cu/Al 2O3 resistive memory switching. Nano Lett.
15(6), 3983 (2015)
24. P.A.M. Dirac, Quantum mechanics of many-electron systems, in
Proceedings of the Royal Society of London A: Mathematical,
Physical and Engineering Sciences, V ol. 123 (The Royal Society,
1929), p. 714
25. J.J. de Pablo, W.A. Curtin, Multiscale modeling in advanced
materials research: Challenges, novel methods, and emerging
applications. MRS Bull. 32(11), 905 (2007)
26. K. Lejaeghere, G. Bihlmayer, T. Bj ¨o r k m a n ,P .B l a h a ,S .B l¨ugel,
V . Blum, D. Caliste, I.E. Castelli, S.J. Clark, A. Dal Corso et
al, Reproducibility in density functional theory calculations of
solids. Science 351(6280), 3000 (2016)
27. A. Shekhar, K. Nomura, R.K. Kalia, A. Nakano, P. Vashishta,
Nanobubble collapse on a silica surface in water: Billion-
atom reactive molecular dynamics simulations. Phys. Rev. Lett.
111(18), 184503 (2013)
28. A.M. Cuiti ˜no, L. Stainier, G. Wang, A. Strachan, T. ˙Ca˙gin,
W.A. Goddard, M. Ortiz, A multiscale approach for modeling
crystalline solids. J. Computer-Aided Mater. Des. 8(2-3), 127
(2001)
29. M. Koslowski, A. Strachan, Uncertainty propagation in a multi-
scale model of nanocrystalline plasticity. Reliab. Eng. Syst. Saf.
96(9), 1161 (2011)
30. M. Ortiz, A.M. Cuitino, J. Knap, M. Koslowski, Mixed
atomistic–continuum models of material behavior The art of tran-
scending atomistics and informing continua. MRS Bull. 26(03),
216 (2001)
31. W.A. Curtin, R.E. Miller, Atomistic/continuum coupling in com-
putational materials science. Model. Simul. Mater. Sci. Eng.
11(3), R33 (2003)
32. J. Guo, S. Datta, M. Lundstrom, M.P. Anantam, Toward multi-
scale modeling of carbon nanotube transistors. Int. J. Multiscale
Comput. Eng. 2(2) (2004)
33. R.P. Vedula, S. Palit, M.A. Alam, A. Strachan, Role of atomic
variability in dielectric charging: A first-principles-based multi-
scale modeling study. Phys. Rev. B 88(20), 205204 (2013)
34. R.A. Austin, N.R. Barton, J.E. Reaugh, L.E. Fried, Direct numer-
ical simulation of shear localization and decomposition reactions
in shock-loaded HMX crystal. J. Appl. Phys. 117(18), 185902
(2015)
35. L. Goux, P. Czarnecki, Y .Y . Chen, L. Pantisano, X.P. Wang, R.
Degraeve, B. Govoreanu, M. Jurczak, D.J. Wouters, L. Altimime,
Evidences of oxygen-mediated resistive-switching mechanism in
TiN\ HfO
2\ Pt cells. Appl. Phys. Lett. 97(24), 243509 (2010)
36. K. Seo, I. Kim, S. Jung, M. Jo, S. Park, J. Park, J. Shin,
K.P. Biju, J. Kong, K. Lee et al, Analog memory and spike-
timing-dependent plasticity characteristics of a nanoscale tita-
nium oxide bilayer resistive switching device. Nanotechnology
22(25), 254023 (2011)
37. J. Song, D. Lee, J. Woo, Y . Koo, E. Cha, S. Lee, J. Park, K. Moon,
S.H. Misha, A. Prakash, et al., Effects of reset current over-
shoot and resistance state on reliability of RRAM. IEEE Electron
Device Lett. 35(6), 636 (2014)
38. O. Kavehei, E. Linn, L. Nielen, S. Tappertzhofen, E. Skafidas,
I. Valov, R. Waser, An associative capacitive network based
on nanoscale complementary resistive switches for memory-
intensive computing. Nanoscale 5(11), 5119 (2013)
39. S. Yu, H.S.P. Wong, Characterization and modeling of the con-
duction and switching mechanisms of HfOx based RRAM, in
MRS Proceedings, V ol. 1631 (Cambridge Univ Press, Cam-
bridge, 2014)
40. L. Zhang, A. Redolfi, C. Adelmann, S. Clima, I.P. Radu, Y .Y .
Chen, D.J. Wouters, G. Groeseneken, M. Jurczak, B. Govore-
anu, Ultrathin metal/amorphous-silicon/metal diode for bipolar
RRAM selector applications. IEEE Electron Device Lett. 35(2),
199 (2014)
41. B. Magyari-K ¨ope, Y . Nishi, Resistive memories. Intelligent Inte-
grated Systems: Devices, Technologies, and Architectures, ed. by
S. Deleonibus, Pan Stanford Series on Intelligent Nanosystems,
vol. 1, p. 325 (2014)
42. S.G. Park, B. Magyari-K ¨ope, Y . Nishi, Electronic correlation
effects in reduced rutile TiO
2 within the LDA+U method. Phys.
Rev. B 82(11), 115109 (2010)
43. S.G. Park, B. Magyari-K ¨ope, Y . Nishi, Theoretical study of the
resistance switching mechanism in rutile TiO 2−x for ReRAM:
the role of oxygen vacancies and hydrogen impurities, in
2011 IEEE Symposium on VLSI Technology-Digest of Technical
Papers, (2011)
44. L. Zhao, S.G. Park, B. Magyari-K ¨ope, Y . Nishi, Dopant selection
rules for desired electronic structure and vacancy formation char-
acteristics of TiO
2 resistive memory. Appl. Phys. Lett. 102(8),
083506 (2013)
45. L. Zhao, S.G. Park, B. Magyari-K ¨ope, Y . Nishi, First principles
modeling of charged oxygen vacancy filaments in reduced TiO2–
implications to the operation of non-volatile memory devices.
Math. Comput. Model. 58(1), 275 (2013)
46. K. Kamiya, M.Y . Yang, S.G. Park, B. Magyari-K ¨ope, Y . Nishi,
M. Niwa, K. Shiraishi, ON-OFF switching mechanism of
resistive–random–access–memories based on the formation and
disruption of oxygen vacancy conducting channels. Appl. Phys.
Lett. 100(7), 073502 (2012)
47. K. Kamiya, M.Y . Yang, T. Nagata, S.G. Park, B. Magyari-
K¨ope, T. Chikyow, K. Yamada, M. Niwa, Y . Nishi, K. Shiraishi,
Generalized mechanism of the resistance switching in binary-
oxide-based resistive random-access memories. Phys. Rev. B
87(15), 155201 (2013)
48. D. Duncan, B. Magyari-K ¨ope, Y . Nishi, Hydrogen doping in
HfO
2 resistance change random access memory. Appl. Phys.
Lett. 108(4), 043501 (2016)
49. D. Duncan, B. Magyari-K ¨ope, Y . Nishi, Filament-induced
anisotropic oxygen vacancy diffusion and charge trapping effects
in hafnium oxide RRAM. IEEE Electron Device Lett. 37(4), 400
(2016)
50. S.R. Bradley, A.L. Shluger, G. Bersuker, Electron-injection-
assisted generation of oxygen vacancies in monoclinic HfO
2.
Phys. Rev. Appl. 4(6), 064008 (2015)
51. S.R. Bradley, G. Bersuker, A.L. Shluger, Modelling of oxygen
vacancy aggregates in monoclinic HfO 2: Can they contribute to
conductive filament formation? J. Phys. Condens. Matter27(41),
415401 (2015)
52. P. Hohenberg, W. Kohn, Inhomogeneous electron gas. Phys. Rev.
136(3B), B864 (1964)
53. W. Kohn, L.J. Sham, Self-consistent equations including exchange
and correlation effects. Phys. Rev. 140(4A), A1133 (1965)

---

<!-- PDF page 20; local proof page count -->

J Electroceram
54. S.L. Dudarev, G.A. Botton, S.Y . Savrasov, C.J. Humphreys, A.P.
Sutton, Electron-energy-loss spectra and the structural stability
of nickel oxide: An LSDA+u study. Phys. Rev. B 57(3), 1505
(1998)
55. T.T. Jiang, Q.Q. Sun, Y . Li, J.J. Guo, P. Zhou, S.J. Ding, D.W.
Zhang, Towards the accurate electronic structure descriptions of
typical high-constant dielectrics. J. Phys. D: Appl. Phys. 44(18),
185402 (2011)
56. M. Jain, J.R. Chelikowsky, S.G. Louie, Quasiparticle excitations
and charge transition levels of oxygen vacancies in hafnia. Phys.
Rev. Lett. 107(21), 216803 (2011)
57. S.J. Clark, L. Lin, J. Robertson, On the identification of the oxy-
gen vacancy in HfO
2. Microelectron. Eng. 88(7), 1464 (2011)
58. K.H. Xue, L.R.C. Blaise, P. Fonseca, G. Molas, E. Vianello, B.
Traore, B. De Salvo, G. Ghibaudo, Y . Nishi, Grain boundary
composition and conduction in HfO 2: An ab initio study. Appl.
Phys. Lett. 102(20), 201908 (2013)
59. G. Kresse, J. Furthm ¨uller, Efficient iterative schemes for ab initio
total-energy calculations using a plane-wave basis set. Phys. Rev.
B 54(16), 11169 (1996)
60. E.P. Bl ¨ochl, Projector augmented-wave method. Phys. Rev. B
50(24), 17953 (1994)
61. B. Xiao, S. Watanabe, Interface structure in Cu/Ta 2O5/Pt resis-
tance switch A first-principles study. ACS Appl. Mater. Inter-
faces 7(1), 519 (2015)
62. A. O’Hara, G. Bersuker, A.A. Demkov, Assessing hafnium on
hafnia as an oxygen getter. J. Appl. Phys.115(18), 183703 (2014)
63. L. Zhao, S. Clima, B. Magyari-K ¨ope, M. Jurczak, Y . Nishi,
Ab initio modeling of oxygen-vacancy formation in doped-
HfOx RRAM: Effects of oxide phases, stoichiometry, and dopant
concentrations. Appl. Phys. Lett. 107(1), 013504 (2015)
64. D. Ning, P. Hua, W. Wei, Effects of different dopants on switch-
ing behavior of HfO
2-based resistive random access memory.
Chin. Phys. B 23(10), 107306 (2014)
65. Y .S. Chen, B. Chen, B. Gao, L.F. Liu, X.Y . Liu, J.F. Kang,
Well controlled multiple resistive switching states in the Al local
doped HfO
2 resistive random access memory device. J. Appl.
Phys. 113(16), 164507 (2013)
66. L. Goux, J.Y . Kim, B. Magyari-K ¨ope, Y . Nishi, A. Redolfi, M.
Jurczak, H-treatment impact on conductive-filament formation
and stability in Ta2O5-based resistive-switching memory cells. J.
Appl. Phys. 117(12), 124501 (2015)
67. Z. Yuanyang, W. Jiayu, X. Jianbin, Y . Fei, L. Qi, D. Yuehua,
Metal dopants in HfO 2–based RRAM: first principle study. J.
Semicond. 35(4), 042002 (2014)
68. B. Gao, H.W. Zhang, S. Yu, B. Sun, L.F. Liu, X.Y . Liu, Y .
Wang, R.Q. Han, J.F. Kang, B. Yu et al, Oxide-based RRAM:
Uniformity improvement using a new material-oriented method-
ology, in 2009 IEEE Symposium on VLSI Technology-Digest of
Technical Papers, V ol. 30 (IEEE, 2009)
69. C.S. Peng, W.Y . Chang, Y .H. Lee, M.H. Lin, F. Chen, M.J. Tsai,
Improvement of resistive switching stability of HfO
2 films with
Al doping by atomic layer deposition. Electrochem. Solid-State
Lett. 15(4), H88 (2012)
70. Z. Wang, W.G. Zhu, A.Y . Du, L. Wu, Z. Fang, X.A. Tran, W.J.
Liu, K.L. Zhang, H.Y . Yu, Highly uniform, self-compliance, and
forming-free ALD-based RRAM with Ge doping. IEEE Trans.
Electron Devices 59(4), 1203 (2012)
71. D. Panda, C.Y . Huang, T.Y . Tseng, Resistive switching character-
istics of nickel silicide layer embedded HfO
2 film. Appl. Phys.
Lett. 100(11), 112901 (2012)
72. S. Kim, D. Lee, J. Park, S. Jung, W. Lee, J. Shin, J. Woo, G.
Choi, H. Hwang, Defect engineering: reduction effect of hydro-
gen atom impurities in HfO
2-based resistive-switching memory
devices. Nanotechnology 23(32), 325702 (2012)
73. H. Xie, Q. Liu, Y . Li, H. Lv, M. Wang, X. Liu, H. Sun, X. Yang,
S. Long, S. Liu, et al., Nitrogen-induced improvement of resistive
switching uniformity in a HfO
2-based RRAM device. Semicond.
Sci. Technol. 27(12), 125008 (2012)
74. G. Palma, E. Vianello, O. Thomas, H. Oucheikh, S. Onkaraiah,
A. Toffoli, C. Carabasse, G. Molas, B. De Salvo, A novel HfO2-
GeS2-Ag based conductive bridge RAM for reconfigurable logic
applications, in Proceedings of the European Solid-State Device
Research Conference (ESSDERC), (2013), p. 2013
75. K.L. Pey, N. Raghavan, X. Wu, M. Bosman, X.X. Zhang, K. Li,
Spatial correlation of conductive filaments for multiple switch-
ing cycles in CBRAM, in 2014 IEEE International Conference
on Electron Devices and Solid-State Circuits (EDSSC) (IEEE,
2014), p. 1
76. T.L. Tsai, H.Y . Chang, F.S. Jiang, T.Y . Tseng, Impact of post-
oxide deposition annealing on resistive switching in HfO
2–based
oxide RRAM and conductive-bridge RAM devices. IEEE Elec-
tron Device Lett. 36(11), 1146 (2015)
77. M.A. Wood, M.J. Cherukara, E.M. Kober, A. Strachan, Ultra-
fast chemistry under nonequilibrium conditions and the shock
to deflagration transition at the nanoscale. J. Phys. Chem. C
119(38), 22008 (2015)
78. M.A. Wood, A. Strachan, Nonequilibrium reaction kinetics in
molecular solids. J. Phys. Chem. C 120(1), 542 (2015)
79. E.C. Neyts, A.C.T. van Duin, A. Bogaerts, Insights in the plasma-
assisted growth of carbon nanotubes through atomic scale simu-
lations: Effect of electric field. J. Am. Chem. Soc. 134(2), 1256
(2011)
80. T.P. Senftle, S. Hong, M.M. Islam, S.B. Kylasa, Y . Zheng, Y .K.
Shin, C. Junkermeier, R. Engel-Herbert, M.J. Janik, H. Aktulga,
T. Verstraelen, A. Grama, A.C.T. van Duin, The reaxFF reactive
force-field: Development, applications and future directions. Npj
Computational Materials 2, 15011 (2016)
81. S.J. Stuart, A.B. Tutein, J.A. Harrison, A reactive potential for
hydrocarbons with intermolecular interactions. J. Chem. Phys.
112(14), 6472 (2000)
82. S.W. Rick, S.J. Stuart, B.J. Berne, Dynamical fluctuating charge
force fields: Application to liquid water. J. Chem. Phys. 101(7),
6141 (1994)
83. A.C.T. Van Duin, S. Dasgupta, F. Lorant, W.A. Goddard, Reaxff:
A reactive force field for hydrocarbons. J. Phys. Chem. A
105(41), 9396 (2001)
84. F.H. Stillinger, T.A. Weber, Computer simulation of local order
in condensed phases of silicon. Phys. Rev. B 31(8), 5262
(1985)
85. J. Tersoff, Modeling solid-state chemistry: Interatomic potentials
for multicomponent systems. Phys. Rev. B 39(8), 5566 (1989)
86. W.J. Mortier, S.K. Ghosh, S. Shankar, Electronegativity-
equalization method for the calculation of atomic charges in
molecules. J. Am. Chem. Soc. 108(15), 4315 (1986)
87. K.A. Rappe, W.A. Goddard III, Charge equilibration for molec-
ular dynamics simulations. J. Phys. Chem. 95(8), 3358 (1991)
88. A.C. Antony, S.A. Akhade, T. Liang, M.J. Janik, J.K. Maranas,
S.B. Sinnott, Simulating an applied voltage in molecular dynam-
ics using charge optimized many body (COMB3) potentials. ECS
Trans. 69(1), 103 (2015)
89. K.D. Nielson, A.C.T. van Duin, J. Oxgaard, W.Q. Deng, W.A.
Goddard, Development of the ReaxFF reactive force field for
describing transition metal catalyzed reactions, with applica-
tion to the initial stages of the catalytic formation of carbon
nanotubes. J. Phys. Chem. A 109(3), 493 (2005)
90. A.C.T. Van Duin, A. Strachan, S. Stewman, Q. Zhang, X.
Xu, W.A. Goddard, ReaxFFSio reactive force field for silicon
and silicon oxide systems. J. Phys. Chem. A 107(19), 3803
(2003)

---

<!-- PDF page 21; local proof page count -->

J Electroceram
91. K. Chenoweth, A.C.T. Van Duin, W.A. Goddard, ReaxFF reac-
tive force field for molecular dynamics simulations of hydrocar-
bon oxidation. J. Phys. Chem. A 112(5), 1040 (2008)
92. N. Onofrio, D. Guzman, A. Strachan, The dynamics of copper
intercalated molybdenum ditelluride. J. Chem. Phys. 145(19),
194702 (2016)
93. N. Onofrio, D. Guzman, A. Strachan, Atomic origin of ultrafast
resistance switching in nanoscale electrometallization cells. Nat.
Mater. 14(4), 440 (2015)
94. N. Onofrio, A. Strachan, V oltage equilibration for reactive atom-
istic simulations of electrochemical processes. J. Chem. Phys.
143(5), 054109 (2015)
95. N. Onofrio, A. Strachan, Alejandro Strachan research group.
https://nanohub.org/groups/strachangroup/lammpsmodules
(2015)
96. D. Bedrov, G.D. Smith, A.C.T. van Duin, Reactions of singly-
reduced ethylene carbonate in lithium battery electrolytes: a
molecular dynamics simulation study using the ReaxFF. J. Phys.
Chem. A 116(11), 2978 (2012)
97. J.T. Su, W.A. Goddard III, Excited electron dynamics modeling
of warm dense matter. Phys. Rev. Lett. 99(18), 185003 (2007)
98. S. Kale, J. Herzfeld, S. Dai, M. Blank, Lewis-inspired represen-
tation of dissociable water in clusters and Grotthuss chains. J.
Biol. Phys. 38(1), 49 (2012)
99. R.A. Nistor, J.G. Polihronov, M.H. M ¨user, N.J. Mosey, A gen-
eralization of the charge equilibration method for nonmetallic
materials. J. Chem. Phys. 125(9), 094108 (2006)
100. M.T. Knippenberg, P.T. Mikulski, K.E. Ryan, S.J. Stuart, G. Gao,
J.A. Harrison, Bond-order potentials with split-charge equilibra-
tion: Application to C-, H-, and O-containing systems. J. Chem.
Phys. 136(16), 164701 (2012)
101. W.B. Dapp, M.H. M ¨user, Redox reactions with empirical poten-
tials: Atomistic battery discharge simulations. J. Chem. Phys.
139(6), 064106 (2013)
102. M.M. Islam, G. Kolesov, T. Verstraelen, E. Kaxiras, A.C.T. van
Duin, eReaxFF: A pseudoclassical treatment of explicit electrons
within reactive force field simulations. J. Chem. Theory Comput.
12(8), 3463 (2016)
103. J.T. Su, W.A. Goddard III, The dynamics of highly excited elec-
tronic systems: Applications of the electron force field. J. Chem.
Phys. 131(24), 244501 (2009)
104. T. Verstraelen, P.W. Ayers, V . Van Speybroeck, M. Waroquier,
ACKS2: Atom-condensed Kohn-Sham DFT approximated to
second order. J. Chem. Phys. 138(7), 074108 (2013)
105. J.C. Rienstra-Kiracofe, G.S. Tschumper, H.F. Schaefer, S. Nandi,
G.B. Ellison, Atomic and molecular electron affinities: photo-
electron experiments and theoretical computations. Chem. Rev.
102(1), 231 (2002)
106. K.D. Jordan, P.D. Burrow, Studies of the temporary anion states
of unsaturated hydrocarbons by electron transmission spec-
troscopy. Acc. Chem. Res. 11(9), 341 (1978)
107. M.M. Islam, A.C.T. van Duin, Reductive decomposition reac-
tions of ethylene carbonate via explicit electron transfer from
lithium An eReaxFF molecular dynamics study. J. Phys. Chem.
C 120(48), 27128 (2016)
108. S.P. Kim, A.C.T. Van Duin, V .B. Shenoy, Effect of electrolytes
on the structure and evolution of the solid electrolyte interphase
(SEI) in Li-ion batteries: A molecular dynamics study. J. Power.
Sources 196(20), 8590 (2011)
109. Y . Bernard, V .T. Renard, P. Gonon, V . Jousseaume, Back-end-
of-line compatible conductive bridging RAM based on Cu and
SiO
2. Microelectron. Eng. 88(5), 814 (2011)
110. N. Onofrio, D. Guzman, A. Strachan, Atomistic simulations
of electrochemical metallization cells: Mechanisms of ultra-fast
resistance switching in nanoscale devices Nanoscale (2016)
111. I. Valov, G. Staikov, Nucleation and growth phenomena in nano-
sized electrochemical systems for resistive switching memories.
J. Solid State Electrochem. 17(2), 365 (2013)
112. A. Belmonte, U. Celano, A. Redolfi, A. Fantini, R. Muller, W.
Vandervorst, M. Houssa, M. Jurczak, L. Goux, Analysis of the
excellent memory disturb characteristics of a Hourglass-shaped
filament in Al
2O3/Cu-based CBRAM devices. IEEE Trans. Elec-
tron Devices 62(6), 2007 (2015)
113. D.M. Guzman, N. Onofrio, A. Strachan, Stability and migration
of small copper clusters in amorphous dielectrics. J. Appl. Phys.
117(19), 195702 (2015)
114. S.C. Pandey, R. Meade, G.S. Sandhu, Cu impurity in insu-
lators and in metal-insulator-metal structures: Implications for
resistance-switching random access memories. J. Appl. Phys.
117(5), 054504 (2015)
115. Z. Jiang, M. Povolotskyi, N. Onofrio, D. Guzman, D. Lemus, S.
Perez, J. Bermeo, A. Strachan, G. Klimeck, Multi-scale quantum
simulations of conductive bridging RAM. IWCE (2015)
116. S.C. Chae, J.S. Lee, S. Kim, S.B. Lee, S.H. Chang, C. Liu,
B. Kahng, H. Shin, D.W. Kim, C.U. Jung, S. Seo, M.J. Lee,
T.W. Noh, Random circuit breaker network model for unipolar
resistance switching. Adv. Mater. 20(6), 1154 (2008)
117. N. Xu, B. Gao, L.F. Liu, B. Sun, X.Y . Liu, R.Q. Han, J.F. Kang,
B. Yu, A unified physical model of switching behavior in oxide-
based RRAM, in 2008 IEEE Symposium on VLSI Technology-
Digest of Technical Papers, (2008), pp. 100–101
118. X. Guan, S. Yu, H.S.P. Wong, On the switching parameter vari-
ation of metal-oxide RRAM – Part I: Physical modeling and
simulation methodology. IEEE Trans. Electron Devices 59(4),
1172 (2012)
119. L. Larcher, A. Padovani, O. Pirrotta, L. Vandelli, G. Bersuker,
Microscopic understanding and modeling of HfO
2 RRAM
device physics, in 2012 IEEE International Electron Devices
Meeting (IEDM) Technical Digest, (2012), p. 474
120. S. Qin, Z. Liu, G. Zhang, J. Zhang, Y . Sun, H. Wu, H. Qian, Z.
Yu, Atomistic study of dynamics for metallic filament growth in
conductive-bridge random access memory. Phys. Chem. Chem.
Phys. 17, 8627 (2015)
121. J. Guy, G. Molas, P. Blaise, M. Bernard, A. Roule, G. Le
Carval, V . Delaye, A. Toffoli, G. Ghibaudo, F. Clermidy, B.
De Salvo, L. Perniola, Investigation of forming, SET, and data
retention of conductive-bridge random-access memory for stack
optimization. IEEE Trans. Electron Devices 62, 3482 (2015)
122. S. Menzel, P. Kaupmann, R. Waser, Understanding filamen-
tary growth in electrochemical metallization memory cells
using Kinetic Monte Carlo simulations. Nanoscale 7, 12673
(2015)
123. F. Pan, S. Yin, V . Subramanian, A detailed study of the forming
stage of an electrochemical resistive switching memory by KMC
simulation. IEEE Electron Device Lett. 32, 949 (2011)
124. A. Padovani, L. Larcher, O. Pirrotta, L. Vandelli, G. Bersuker,
Microscopic modeling of HfOx RRAM operations From forming
to switching. IEEE Trans. Electron Devices 62, 1998 (2015)
125. J. Hur, D. Lee, S. Jeon, A Monte Carlo simulation for bipo-
lar resistive memory switching in large band-gap oxides. Appl.
Phys. Lett. 203504, 107 (2015)
126. D. Ielmini, Modeling the universal set/reset characteristics of
bipolar RRAM by field-and temperature-driven filament growth.
IEEE Trans. Electron Devices 58(12), 4309 (2011)
127. S. Yu, H.S.P. Wong, Compact modeling of conducting-bridge
random-access memory (CBRAM). IEEE Trans. Electron
devices 58(5), 1352 (2011)
128. S. Lin, L. Zhao, J. Zhang, H. Wu, Y . Wang, H. Qian, Z.
Yu, Electrochemical simulation of filament growth and disso-
lution in conductive-bridging RAM (CBRAM) with cylindrical

---

<!-- PDF page 22; local proof page count -->

J Electroceram
coordinates, in 2012 IEEE International Electron Devices Meet-
ing (IEDM) Technical Digest, (2012), p. 593
129. S. Ambrogio, S. Balatti, D.C. Gilmer, D. Ielmini, Analytical
modeling of oxide-based bipolar resistive memories and com-
plementary resistive switches. IEEE Trans. Electron Devices 61,
2378–2386 (2014)
130. S. Yu, Y . Wu, H.S.P. Wong, Investigating the switching dynam-
ics and multilevel capability of bipolar metal oxide resistive
switching memory. Appl. Phys. Lett. 98(10), 103514 (2011)
131. Z. Fang, H.Y . Yu, W.J. Fan, G. Ghibaudo, J. Buckley, B. DeSalvo,
X. Li, X.P. Wang, G.Q. Lo, D.L. Kwong, Current conduction
model for oxide-based resistive random access memory verified
by low-frequency noise analysis. IEEE Trans. Electron Devices
60(3), 1272 (2013)
132. S. Lv, H. Wang, J. Zhang, J. Liu, L. Sun, Z. Yu, An analytical
model for the forming process of conductive-bridge resistive-
switching random-access memory. IEEE Trans. Electron Devices
61, 3166–3171 (2014)
133. S. Menzel, S. Tappertzhofen, R. Waser, I. Valov, Switching kinet-
ics of electrochemical metallization memory cells. Phys. Chem.
Chem. Phys. 15, 6945 (2013)
134. P. Sun, L. Li, N. Lu, Y . Li, M. Wang, H. Xie, S. Liu, M. Liu, Phys-
ical model of dynamic Joule heating effect for reset process in
conductive-bridge random access memory. J. Comput. Electron.
13, 432 (2014)
135. R. Degraeve, A. Fantini, S. Clima, B. Govoreanu, L. Goux,
Y .Y . Chen, D.J. Wouters, P. Roussel, G.S. Kar, G. Pourtois, S.
Cosemans, J.A. Kittl, G. Groeseneken, M. Jurczak, L. Altim-
ime, Dynamic ’hour glass’ model for SET and RESET in HfO
2
RRAM, in 2012 IEEE Symposium on VLSI Technology-Digest of
Technical Papers(IEEE, 2012), p. 75
136. S. Ambrogio, S. Balatti, A. Cubeta, A. Calderoni, N. Ramaswamy,
D. Ielmini, Statistical fluctuations in HfOx resistive-switching
m e m o r yP a r tI-S e t / R e s e tv a r i a bility. IEEE Trans. Electron
Devices 61(8), 2912 (2014)
137. N. Raghavan, R. Degraeve, A. Fantini, L. Goux, S. Strangio, B.
Govoreanu, D.J. Wouters, G. Groeseneken, M. Jurczak, Micro-
scopic origin of random telegraph noise fluctuations in aggres-
sively scaled RRAM and its impact on read disturb variability, in
2013 IEEE International Reliability Physics Symposium (IRPS)
(IEEE, 2013), pp. 5E–3
138. A. Fantini, L. Goux, R. Degraeve, D.J. Wouters, N. Raghavan, G.
Kar, A. Belmonte, Y .-Y . Chen, B. Govoreanu, M. Jurczak, Intrin-
sic switching variability in HfO2 RRAM, in Memory Workshop
(IMW), 2013 5th IEEE International (IEEE, 2013), p. 30
139. I.G. Baek, M.S. Lee, S. Seo, M.J. Lee, D.H. Seo, D.S. Suh,
J.C. Park, S.O. Park, H.S. Kim, I.K. Yoo, et al., Highly scalable
nonvolatile resistive memory using simple binary oxide driven
by asymmetric unipolar voltage pulses, in 2004 IEEE Inter-
national Electron Devices Meeting (IEDM) Technical Digest
(IEEE, 2004), p. 587
140. I.G. Baek, D.C. Kim, M.J. Lee, H.J. Kim, E.K. Yim, M.S. Lee,
J.E. Lee, S.E. Ahn, S. Seo, J.H. Lee, et al., Multi-layer cross-
point binary oxide resistive memory (oxrram) for post-NAND
storage application, in 2005 IEEE International Electron Devices
Meeting (IEDM) Technical Digest (IEEE, 2005), p. 750
141. G. Indiveri, B. Linares-Barranco, R. Legenstein, G. Deligeor-
gis, T. Prodromakis, Integration of nanoscale memristor synapses
in neuromorphic computing architectures. Nanotechnology 24,
384010 (2013)
142. M. Suri, D. Querlioz, O. Bichler, G. Palma, E. Vianello, D.
Vuillaume, C. Gamrat, B. DeSalvo, Bio-inspired stochastic com-
puting using binary CBRAM synapses. IEEE Trans. Electron
Devices 60(7), 2402 (2013)
143. B.L. Jackson, B. Rajendran, G.S. Corrado, M. Breitwisch, G.W.
Burr, R. Cheek, K. Gopalakrishnan, S. Raoux, C. Rettner, A.
Padilla, A.G. Schrott, R.S. Shenoy, B.N. Kurdi, C.H. Lam,
D.S. Modha, Nanoscale electronic synapses using phase change
devices. ACM J. Emerg. Technol. Comput. Syst.9(2), 1 (2013)
144. M. Youssef, B. Yildiz, Predicting self diffusion in metal oxides
from first principles: The case of oxygen in tetragonal ZrO
2.
Phys. Rev. B 89(2), 024105 (2014)
145. M. Panzer, M. Shandalov, J. Rowlette, Y . Oshima, Y .W. Chen, P.
McIntyre, K. Goodson, Thermal properties of ultrathin hafnium
oxide gate dielectric films. IEEE Electron Device Lett. 30(12),
1269 (2009)
146. S. Larentis, F. Nardi, S. Balatti, D.C. Gilmer, D. Ielmini, Resis-
tive switching by voltage-driven ion migration in bipolar RRAM
– Part II Modeling. IEEE Trans. Electron Devices 59(9), 2468
(2012)
147. S. Kim, S.J. Kim, K.M. Kim, S.R. Lee, M. Chang, E. Cho, Y .B.
Kim, C.J. Kim, U.I. Chung, I.K. Yoo, Physical electro-thermal
model of resistive switching in bi-layered resistance-change
memory. Sci. Rep. 3, 1680 (2013)
148. F. Nardi, D. Ielmini, C. Cagli, S. Spiga, M. Fanciulli, L. Goux,
D.J. Wouters, Control of filament size and reduction of reset
current below 10μA in NiO resistance switching memories.
Solid-State Electron. 58(1), 42 (2011)
149. A. Marchewka, B. Roesgen, K. Skaja, H. Du, C.L. Jia, J. Mayer,
V . Rana, R. Waser, S. Menzel, Nanoionic resistive switching
memories: On the physical nature of the dynamic reset process.
Advanced Electronic Materials 2, 1500233 (2016)
150. J.S. Lee, S.B. Lee, B. Kahng, T.W. Noh, Two opposite hysteresis
curves in semiconductors with mobile dopants. Appl. Phys. Lett.
102, 253503 (2013)
151. S. Ambrogio, S. Balatti, A. Cubeta, A. Calderoni, N.
Ramaswamy, D. Ielmini, Statistical fluctuations in HfOx
resistive-switching memory Part II – random telegraph noise.
IEEE Trans. Electron Devices 61(8), 2920 (2014)
152. S. Ambrogio, S. Balatti, V . McCaffrey, D.C. Wang, D. Ielmini,
Noise-induced resistance broadening in resistive switching mem-
ory – Part II Array statistics. IEEE Trans. Electron Devices
62(11), 3812 (2015)
153. K.S. Ralls, W.J. Skocpol, L.D. Jackel, R.E. Howard, L.A. Fetter,
R.W. Epworth, D.M. Tennant, Discrete resistance switching in
submicrometer silicon inversion layers Individual interface traps
and low-frequency (1 f) noise. Phys. Rev. Lett.52(3), 228 (1984)
154. R. Soni, P. Meuffels, A. Petraru, M. Weides, C. K ¨ugeler, R.
Waser, H. Kohlstedt, Probing Cu doped Ge
0.3Se0.7 based resis-
tance switching memory devices with random telegraph noise. J.
Appl. Phys. 107(2), 024517 (2010)
155. D. Veksler, G. Bersuker, L. Vandelli, A. Padovani, L. Larcher, A.
Muraviev, B. Chakrabarti, E. V ogel, D.C. Gilmer, P.D. Kirsch,
Random telegraph noise (RTN) In scaled RRAM devices, in
2013 IEEE International Reliability Physics Symposium (IRPS)
(IEEE, p. 2013