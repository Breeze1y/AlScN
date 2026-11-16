<!-- PDF page 1 -->

Device Simulation and Modeling of RRAM Switching Processes
Zhiping Yu and Jinyu Zhang
Institute of Microelectronics, Tsinghua University, Beijing 100084, China
yuzhip@tsinghua.edu.cn
Abstract
On the emerging nonvolatile memory (NVM) devices,
resistive-switching RAMs (ReRAM or RRAM), are con-
sidered as promising next-generation structures. All
types of RRAMs (TMO, CBRAM, CMO, and etc.) work
under the same principle of conductive-ﬁlament forming
and rupturing in the switching layer of the two-terminal
devices. This switching mechanism is to diﬀerentiate
R R A Mf r o mo t h e rt y p e so fN V M ss u c ha sP C M( p h a s e -
change memory) and STT-MRAM (spin-transfer torque
magnetic RAM).
This talk describes the device simulation and mod-
eling of switching process in RRAMs and in particular,
the conductive-bridge RAMs (CBRAMs).
Speciﬁcally, 3D kinetic Monte-Carlo (3D-KMC) sim-
ulation and a compact model for the forming of ﬁlament
will be presented. A review of state-of-the-art in this re-
search ﬁeld is provided.
1 Introduction
In view of Intel’s announcement in July 2015 of up-
coming 3D XPoint memory, which they term as ”the
ﬁrst new class of memory since the introduction of NAND
ﬂash 25 years ago” [1], the search for the exact mecha-
nism responsible for switching in NVMs in general and
in RRAMs in particular is on the high notes.
As an example, IBM Z¨u r i c hL a b ,i ni t sNature Com-
munications article [2], measured the crystal growth ve-
locity of a phase change material at a PCM cell and re-
vealed the temperature dependence of the crystal growth
velocity.
Because of the nanoscale and random nature of ﬁl-
aments in RRAM, experimental means are of limited
role, though critical. The physics-based simulation and
modeling on ﬁlament evolution mechanism become es-
sential.
In CBRAMs, the process of ﬁlament growth and rup-
ture) is basically an electrochemical process, but occurs
in a stochastic way. In addition, as in PCM cells, the
temperature plays a critical role in this process.
A cell of two-terminal CBRAM with forming and
dissolving processes is shown schematically in Fig. 1.
We provide two modeling approaches in this review
talk: numerical simulation and in particular kinetic Monte-
Carlo (KMC) in 3D is ﬁrst introduced. Then, a compact
model for the forming process in CBRAM is presented.
Figure 1: Schematic view of dynamic process for CF
forming and dissolving within a CBRAM cell
2 3D Kinetic Monte-Carlo Simulation
The Kinetic Monte-Carlo (KMC) method has suc-
cessfully been employed to model RRAMs [3]. The na-
ture of RRAM is characterized by randomness. For
CBRAMs, the device operation involves chemical reac-
tions and physical migration of ions. KMC is an eﬀective
approach in that it can eﬀectively capture the random
nature of ﬁlament dynamics and in the meantime the
computation time is manageable.
A new 3D KMC program has been developed in au-
thor’s group [4] to capture the four processes in ﬁlament
evolution: forming, thickening, rupturing, and dissolv-
ing. In particular, the dissolution of the conductive ﬁl-
ament (CF) is modeled and simulated with the aid of a
resistor network to incorporate the eﬀect of Joule heat-
ing. The validity of the simulation tool is conﬁrmed by
the simulated reverse proportional relationship between
the on-state resistance and the compliance current as
well as the exponential dependence of forming time over
forming voltage, both are observed in experiments.
3 Compact Modeling of CBRAM
Filament Forming Process
An analytical model for the forming process in CBRAM
is also developed in [5]. The measurable forming time
can be calculated using this model given the bias con-
dition and is veriﬁed through comparison with the ex-
perimental data. The forming time has been shown to
have multiple-slopes in exponential dependence on the
applied voltage, in agreement with measurement. The
model is based on the identiﬁcation of three steps in the
forming process: 1) metal oxidation at the anode; 2)
cation migration toward the cathode; and 3) backward
1

---

<!-- PDF page 2 -->

ﬁlament growth due to electroplating from the cathode.
The accuracy of the model has also been established via
numerical simulation.
The mathematical model for the forming process is
partitioned into ﬁve steps, with schematic shown in Fig. 2.
Closed form is obtained for calculating the forming time.
Figure 2: Schematic view of diﬀerent stages for the CF
growth in the forming process of a CBRAM cell. (a)
Initial migration stage. (b) LFG growth stage. (c) HFG
growth stage. (d) CF reached anode. (e) Lateral growth
stage.
The model can be used in circuit simulation.
4 Other Vacancy-based Models for
TMO-type of RRAMs
Even though all RRAMs have similar ﬁlamental switch-
ing mechanism, TMO (transition-metal)-RRAM, such
as HfO x-based RRAM, has diﬀerent physical process
from CBRAM. The switching characteristics is strongly
correlated with the geometry of CF in TMO-RRAM as a
direct result of generation and recombination of oxygen
vacancies (V
O) in the switching oxide layer.
A physical based analytical model of metal-oxide based
RRAM cell under DC and pulse operation modes has
been developed in Peking University [6]. In this model,
the transport behaviors of oxygen vacancies and oxy-
gen ions, metal conductivity, electron hopping and heat
conduction and the parasitic capacitance and resistance
eﬀects are covered. The developed analytical model is
veriﬁed and calibrated by measured data.
5 Summary
To realize the next-generation of NVM beyond NAND-
ﬂash with cross-bar architecture, both materials and
switching mechanisms for the memory cells, need in-
novation and thorough understanding. Simulation and
modeling are two essential tools in helping optimize the
device structure and memory circuit design. Further-
more, the Monte-Carlo method together with the ther-
mal analysis may be well positioned to resolve the device
variation and reliability issues.
REFERENCES
[1] Peter Clarke, “Intel, Micron launch “bulk-
switching” ReRAM,”EE Times, July 28, 2015.
[2] Abu Sebastian, Manuel Le Gallo, & Daniel Krebs,
“Crystal growth within a phase change memory
cell,”, Nature Communications, p. 1, July 7, 2014
(DOI: 10.1038/ncomms5314)
[3] F. Pan, “Experimental and Simulation Study of Re-
sistive Switches for Memory Applications,” Ph. D.
thesis, UC Berkeley, 2012.
[4] Shi Dong, et al. “Three-Dimensional Kinetic
Monte-Carlo (KMC) Simulation of Conductive-
Bridge Random Access Memory (CBRAM),” pri-
vate communication, Tsinghua, 2016
[5] Shaoli Lv, He Wang, Jinyu Zhang, Jun Liu,
Lingling Sun, and Zhiping Yu, “An Analytical
Model for the Forming Process of Conductive-
Bridge Resistive-Switching Random-Access Mem-
ory (CBRAM),” IEEE T-ED , Vol. 61, No. 9, p.
3166, Sept. 2014.
[6] P. Huang, X. Liu, et al., “A physical based analytic
model of RRAM operation for circuit simulation,”
IEDM, p. 605, 2012.
2