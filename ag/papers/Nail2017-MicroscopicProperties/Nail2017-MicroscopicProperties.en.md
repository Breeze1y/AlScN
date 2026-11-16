<!-- PDF page 1 -->

This article has been accepted for inclusion in a future issue of this journal. Content is final as presented, with the exception of pagination.
IEEE TRANSACTIONS ON ELECTRON DEVICES 1
A Link Between CBRAM Performances and
Material Microscopic Properties Based on
Electrical Characterization and
Atomistic Simulations
Cécile Nail, Gabriel Molas, Philippe Blaise, Benoit Sklenard, Remy Berthier, Mathieu Bernard,
Luca Perniola, Gérard Ghibaudo, and Christophe Vallée,Life Fellow, IEEE
Abstract — In this paper, we investigate the link between
various resistive memory (RRAM) electrical characteristics:
endurance, window margin (WM), and retention. For this
purpose, several RRAMs are characterized using various
resistive layers and bottom electrodes. By focusing on
one technology and optimizing programming conditions
(current, voltage, and time), we establish a tradeoff between
endurance and WM. Then, by changing memory stack,
we demonstrate the correlation between endurance plus
window marging improvement and retention degradation.
Studying this last feature from a material point of view,
we analyze different oxides by density functional theory.
We realize a systematic review for possible exchanges of
species between resistive layer and Cu-based top elec-
trode and study their diffusion. This provides insights
on conductive ﬁlament composition in different stacks.
Combining previous experiments and simulations, we pro-
pose a link between memory characteristics and material
microscopic parameters, through the ion energy migra-
tion barrier. Finally, we extract how endurance, WM, and
retention are correlated to material properties and electrical
parameters in order to choose the suitable material for a
deﬁned application using the RRAM technology.
Index Terms — Ab initio calculation, defect diffusion,
endurance, RRAM tradeoff, thermal stability, window
margin (WM).
I. I NTRODUCTION
R
ESISTIVE memories present promising performances,
such as high density, high speed, and integration in the
BEOL of a logic process. Recently, good characteristics were
reported, making RRAM potential candidates to replace Flash
for both stand alone, storage class applications [1] and embed-
ded products [2], [3]. Endurance over 108 cycles is targeted for
Manuscript received June 23, 2017; revised August 10, 2017; accepted
August 30, 2017. The review of this paper was arranged by Editor J. Kang.
(Corresponding author: Cécile Nail.)
C. Nail is with CEA, LETI, MINATEC Campus, 38054 Grenoble, France
(e-mail: cecile.nail@gmail.com).
G. Molas, P . Blaise, B. Sklenard, R. Berthier, M. Bernard, and
L. Perniola are with CEA, LETI, MINATEC Campus, 38054 Grenoble,
France.
G. Ghibaudo is with IMEP LAHC CNRS, 38016 Grenoble, France.
C. Vallée is with LTM CNRS, 38054 Grenoble, France.
Color versions of one or more of the ﬁgures in this paper are available
online at http://ieeexplore.ieee.org.
Digital Object Identiﬁer 10.1109/TED.2017.2750910
Fig. 1. WM (mean value) as a function of RRAM maximum achieved
cycles reported in the literature depending on maximum temperature
stability for retention after 24 h.
storage class compared with 10 5–106 cycles for stand-alone
applications [4] while for automotive applications tempera-
ture stability over 175 °C is required [5]. In [6], endurance
up to 10 12 cycles was reported with the TiN/GdO/Hf/TiN
stack with a retention below 150 °C. Excellent stability was
demonstrated in [7], where a retention up to 350 °C with
an endurance up to 10
3 cycles was shown. However, it is
very challenging to combine all features, such as having
good cycling, stable retention at high temperature as well as
high enough window margin (WM) depending on the targeted
density. Indeed, intrinsic RRAM variability is a challenge to
avoid WM closure [8]. In [1], a 10-Mb array integration is
demonstrated with a median memory window of three orders
of read current, which leads to an overlap at the −4.6σ
level ( >99.9% of the devices). Fig. 1 shows memory per-
formances data based on experimental results reported in the
literature. It illustrates a general trend where high WM corre-
sponds to low endurance. Meanwhile, there are some outliers
showing more promising characteristics combining high WM
and endurance (blue references). However, these technologies
exhibit lower temperature stability for 24-h retention. In order
to gather a maximum of data, a temperature stability criterion
is taken as the maximum demonstrated temperature with no
resistance variation during 24 h. Based on this observation,
in this paper, we propose to investigate the interrelationship
between these three memory features.
0018-9383 © 2017 IEEE. Personal use is permitted, but republication/redistribution requires IEEE permission.
See http://www.ieee.org/publications_standards/publications/rights/index.html for more information.

---

<!-- PDF page 2 -->

This article has been accepted for inclusion in a future issue of this journal. Content is final as presented, with the exception of pagination.
2 IEEE TRANSACTIONS ON ELECTRON DEVICES
Fig. 2. Schematic of the RRAM families studied in this paper with TEM
cross sections of W/HfO2/CuTex-based CBRAM and Energy-dispersive
X-ray spectroscopy (Hf, W, and Cu) images.
This paper focuses on oxide-based conductive bridge ran-
dom access memory (CBRAM), using Cu-rich-based c onduc-
tive ﬁlament inside an oxide insulating layer [9], sandwiched
between a top electrode (TE) and a bottom electrode (BE)
[10], [11]. As it was previously demonstrated, oxygen vacan-
cies (VO ) assist copper (Cu) diffusion in the RRAM switching
mechanism [12]–[15], and VO and Cu are core to our density
functional theory (DFT) microscopic study.
In this paper, 1T1R oxide-based CBRAM devices are
electrically characterized to hi ghlight the tradeoff between
endurance and WM in one RRAM technology. Then, com-
parison between different stacks is made to illustrate the
correlation between endurance, WM, and retention. To go
further, we investigate the role of copper and oxygen vacancies
in different technologies with the DFT. Focus is put on
HfO
2/CuTex -based CBRAM. Finally, ﬁrst principle calcula-
tions are used to give insights into ﬁlament properties and
composition. We bring out a link between oxide properties and
memory features particularly WM, endurance, and retention.
This allows to build an optim ized memory stack reaching
speciﬁcations of a given memory application.
II. S AMPLE DESCRIPTION
Our measurements are done on different RRAM samples
with a BE/Oxide/CuTe x Gey structure (see Fig. 2). The mem-
ory is integrated in 1T1R conﬁgurations using a CMOS tran-
sistor as a selector. The 3–3.5 nm of amorphous a-Al
2O3 [16]
or 5 nm of amorphous a-HfO 2 [15] or 4 nm of GdO x [17]
was deposited by PVD or ALD on different PVD BEs
(TiN and W). The 20 nm of CuTe
x Gey deposited by PVD
was used as TE. Fig. 2 shows the TEM cross section of
W/HfO2/CuTex RRAM developed in this paper using a plug
structure.
III. E XPERIMENTAL RESULTS
Focusing on TiN/Al 2O3/CuTex samples, we ﬁrst illustrate
the tradeoff between WM and SET programming energy:
increasing the SET current to 100 μA improves the WM.
Beyond that, the WM is degraded due to high electric
stress applied to the dielectric [Fig. 3(a) ], reducing the
OFF-state resistance after RESET operation. Moreover, reduc-
ing the pulse length for the SET and RESET conditions
[see Fig. 3(b) ] reduces the WM, but increases speed and
decreases consumption.
To further study the tradeoff between WM, endurance,
and RESET conditions, we extract WM and the maximum
number of cycles a cell can reach for different RESET
voltages ( V
RESET) setting time pulse at 1 μs. The maximum
number of cycles ( Ncmax) is achieved when the RESET opera-
tion becomes impossible, meaning that the oxide layer was led
to an irreversible breakdown [see Fig. 3(c) and (d)], leading
Fig. 3. (a) Summary of WM behavior depending on SET current Ic and
programming time tprog.L o w Ic offers low WM, while high Ic (150μA)
degrades ROFF. (b) Long tprog improves the WM at the expense of lower
speed and higher consumption. (c) and (d) Endurance characteristics
of TiN/Al 2O3/CuTex RRAM at 100 μA with cycle-to-cycle resistance
distribution (tSET = tRESET = 1μs). Error bars correspond to standard
deviation measured on 100 subsequent cycles.
Fig. 4. a) Typical endurance curves with different VRESET voltages to
extract (b) and (c). Impact of VRESET on TiN/Al2O3/CuTex sample on
(b) ROFF, error bars corresponds to standard deviation on ten cells and
(c) endurance (maximum number of achieved cycles).
to WM closure. Fig. 4(a) and (b) shows VRESET inﬂuence on
ROFF and endurance, respectively. ROFF increases with VRESET
until a threshold value is reached, where the WM is then
degraded. The WM degradation after 2 V can be explained by
a negative SET following the RESET due to defect injection
from the BE [18]. On the other hand, the maximum number of
cycles decreases as V
RESET increases. These ﬁrst results reveal
a tradeoff between WM and endurance: staying out of an early
degradation with VRESET below a threshold value, increasing
WM leads to a degraded endurance.
In order to compare different materials, this paper is
extended to four RRAM technologies. We chose to ﬁx the
programming current to 100 μA, which is close to the opti-
mum programming current, insuring large WM and the high
number of cycles. By focusing on the RESET parameter giving
the best endurance for each technology, we compared their

---

<!-- PDF page 3 -->

This article has been accepted for inclusion in a future issue of this journal. Content is final as presented, with the exception of pagination.
NAIL et al.: LINK BETWEEN CBRAM PERFORMANCES AND MATERIAL MICROSCOPIC PROPERTIES 3
Fig. 5. Endurance and retention performances for (a) W/HfO2/CuTex,
(b) TiN/Al2O3/CuTex, and (c) TiN/GdOx/CuTex memory technologies.
Fig. 6. (a) WM as a function of maximum endurance. For each RRAM
technology, various WMs were achieved changing the RESET conditions.
G r a yl i n ei su s e dt oe x t r a c t(b) for a constant WM. b) Maximum stable
retention temperature after 24-h baking time as a function of maximum
endurance for a given WM (
ROFF/RON ∼ 400).
endurance ( Ncmax) and their best retention (highest stability
temperature T◦
max for 24-h baking time). Results are shown
in Fig. 5 for three of them. W/HfO 2/CuTex shows up to
6.109 cycles with a constant WM while its maximum temper-
ature stability is around 100 °C. TiN/Al 2O3/CuTex shows an
endurance up to 10 7 and a temperature stability up to 200 °C
while TiN/GdOx /CuTex shows an endurance up to 3.103 and a
stable retention up to 300 °C. These three RRAMs are a good
example to illustrate the tradeoff between endurance and reten-
tion: good cycling goes with poor retention (W/HfO 2/CuTex
technology), while the stable r etention is reached at the cost
of limited endurance (TiN/GdOx /CuTex technology). Fig. 6(a)
represents endurance reached by each RRAM technology for
different RESET voltage conditions, implying different WMs.
Focusing on one material, the WM can be increased at the
expense of a degraded endurance. Having compared different
materials (for a given WM), it can be concluded that endurance
is improved as retention is degraded [ Fig. 6(b) ]. Our second
general conclusion is that a tradeoff exists between endurance,
WM, and retention in RRAM technologies.
IV. S
IMULATIONS FRAMEWORK
Atomistic simulations were realized with Siesta code [19]
within DFT [20], [21] using a generalized gradient approxi-
mation (GGA/Perdew, Burke, and Ernzerhof) and Troullier–
Fig. 7. HfO2 and a-HfO 2 structure with different copper interstitial
positions. O in red, Hf in blue, and Cu in orange. Cu i1 and Cu i2 are
used for migration barrier calculation.
Martins pseudopotentials for Al 3s 23p1,H f6 s 25d2,G d
6s25d14f7,O2 s 22p4,C u 3 d 104s1,T e4 d 105s25p4,a n dG e
3d104s24p2 to account for core electrons [22]. A polarized
double-zeta basis set with an energy shift of 50 meV and a
mesh cutoff between 300 and 400 Rydberg with a conver-
gence criterion for the maximum force acting on an atom to
be smaller than 0.02 eV/ ˙A were used. The nudged elastic
band (NEB) calculation was done with 10–20 images, 1
k point, a spring constant of 0.1 eV/ ˙A, a maximum force
of 0.1 eV/ ˙A, and a climbing image. For hafnia, a 96-atom
supercell consisting of eight monoclinic-HfO
2 unit cells was
calculated. This structure is the closest crystalline form of our
amorphous PVD material with a density of 9 .8g ·cm−3 versus
the measurement value of 10 g·cm−3. Different positions were
tested for each defect, especially for interstitials, up to four
inequivalent sites were compared (see Fig. 7).
The band gap obtained at 3.73 eV is corrected [23] in
G0W0 at 5.2 eV versus 5.5 eV measured by ellipsometry.
For the TE, we used a 96-atom supercell consisting of eight
Cu
2Te3Ge unit cells [24]. This structure is the closest crys-
talline form of our amorphous PVD material with a density
of 6.0 g · cm−3 versus 6.1 g · cm−3 (measured). All results
are based on a band alignment between the bandgap of hafnia
and Cu
2Te3Ge, which was done considering: 1) G0W0 gap
correction with an opening of 1 eV at the conduction band
and −0.5 eV at the valence band and 2) electron afﬁnity of
HfO2 at 2 eV [25]. To complete the study, an amorphous
a-HfO2 model with a density of 9.4 g·cm−3 was implemented.
This model was generated using the melt and quench method:
a crystalline monoclinic-HfO2 with Schottky defects was used
for the initial structure to lower density with the eight HfO 2
unit cells considered previously. It was melted at 2800 K for
10 ps and cooled to 300 K for 10 ps with a time step of 1 fs.
A 1 k point was used during the molecular dynamics while
2×2×2k-points were used for the defect calculations. The
band gap is obtained at 3 eV free of defect against 3.73 eV
in the crystalline structure.
V. C
ALCULATION RESULTS
In this section, the objective is to analyze the ﬁlament
composition of the various investigated memory stacks and
to extract the microscopic material properties, in partic-
ular the energy required to diffuse ions during ﬁlament
formation. We studied HfO
′
2s point defects in contact with

---

<!-- PDF page 4 -->

This article has been accepted for inclusion in a future issue of this journal. Content is final as presented, with the exception of pagination.
4 IEEE TRANSACTIONS ON ELECTRON DEVICES
Fig. 8. Formation enthalpy calculation of defects in HfO 2/CuTexGey
system. The Fermi level of three types of metal (Pt, W–Cu, and Ti–Al) is
represented (dashed vertical lines).
aC u T ex -based TE. As we consider exchanges between HfO 2
and CuTe2Ge3, calculations are based on the following equa-
tion for an interstitial of copper charged +q( C u+q
i )
/Delta1H = UHfO2:Cu+q
i
− UHfO2 + UCuTe2Ge3:VCu − UCuTe2Ge3
+ q · (EFermi + /Delta1V + EVB) (1)
where /Delta1H the formation enthalpy of the defect, UDefect
the total energy of the defect, UPristine the total energy of
the initial state, ni and μi the number of chemical formula
and the reference chemical potential for the exchange, q the
charges exchanged in the system, EFermi the Fermi level
of the electrode in considered contact, /Delta1V the electrostatic
energy shift [26] and EVB the highest valence band energy
level. Results are presented in Fig. 8 with the lowest charged
defects. Considering the charge of the species and a TE in
contact based on Cu (vertical orange line), the most favorable
interstitials and vacancies we could obtain are, respectively:
copper interstitial (Cu +1
i ) and oxygen vacancy ( V+2
O ).
Tellurium interstitial and hafnium vacancies are not considered
here because of their high formation enthalpy energy.
Based on these two candidates, we studied the case of
copper introduction in defective hafnia. The reference for the
oxide becomes then a V
O -rich HfO 2 (HfO2 with one VO
defect) instead of stoichiometric HfO 2. VO in HfO2 can come
from: 1) HfO 2 process condition in O-poor atmosphere or
2) O exchange with the TE during forming process requir-
ing 4.8 eV . The results are based on the following equation :
/Delta1H = UHfO2:VO+Cu+q
i
− UHfO2:VO + UCuTe2Ge3:VCu
− UCuTe2Ge3 + q · (EFermi + /Delta1V + EVB). (2)
Cu introduction is favorable in interstitial site while if HfO2
is VO rich, Cu introduction is the easiest in the interstitial site
next to VO . Table I summarizes these results.
To complete these results, the same calculations were per-
formed in amorphous-HfO 2 (see structure Fig. 7 ). Up to six
positions were tested for each defect. Table II reports the most
favorable ones. The results are similar to crystalline HfO 2,
allowing to pursue the calculation only with this structure to
optimize calculation time.
TABLE I
FORMATION ENTHALPY ΔH (eV) B ETWEEN THE TE AND HfO2 EITHER
STOICHIOMETRIC OR INCLUDING OXYGEN VACANCIES
TABLE II
FORMATION ENTHALPY (ΔH) IN a-HfO 2 STOICHIOMETRIC
OR VO-RICH CALCULATED IN GGA C ONSIDERING
EXCHANGES WITH A Cu2Te3Ge TE
Fig. 9. First principle calculations to compute the barrier height to
diffuse from Cu i1 to Cu i2 sites (see Fig. 7): Cu+1
i diffusion barrier in
stoichiometric or defective HfO2 systems.
To go further, Cu migration barrier calculations were
performed. While thermodynamic approach only considers the
initial and ﬁnal energy states of the defects, more insights are
necessary to quantify the energetic barrier height to diffuse one
Cu ion from one site to the other. To this aim, energy migration
barriers were calculated using NEB-based approaches. The
presented results correspond to Cu migration from Cu
i1 site
to Cui2 (see Fig. 7).
The diffusion barrier of Cu +
i in HfO 2 structure was found
close to 2.5 eV ( Fig. 9 ). We also calculated Cu migration
barriers in defective hafnia structures with a VO next to the
Cu migration path. In particular, the diffusion barriers of Cu +
i
in HfO2-VO are evaluated at 1.28 eV . Cu can pass through a
CuO state.
This paper is done for three different materials: Al 2O3,
HfO2,a n dG d 2O3. Table III summarizes the main results.
If thermodynamic and kinetic energies are taken into
account, Al 2O3 and Gd 2O3 show (Cui VO)+1 defect as
the most favorable while HfO 2 shows favorable energies
for (Cui VO)+1 and V+2
O .

---

<!-- PDF page 5 -->

This article has been accepted for inclusion in a future issue of this journal. Content is final as presented, with the exception of pagination.
NAIL et al.: LINK BETWEEN CBRAM PERFORMANCES AND MATERIAL MICROSCOPIC PROPERTIES 5
TABLE III
MIGRATION ENERGY BARRIER , Ea, AND FORMATION ENTHALPY ,ΔH,B ETWEEN THE TE AND DIFFERENT OXIDES EITHER
STOICHIOMETRIC OR INCLUDING OXYGEN VACANCIES (CUiVO DEFECT )
Fig. 10. Emigration impact study on (a) endurance with typical curves
to extract endurance at different WMs and (b) retentions. Schema of an
RRAM mechanism facilitating (c) endurance and (d) retention.
VI. D ISCUSSION
Combining electrical charact erization and DFT results,
we put in evidence a clear correlation between memory
characteristics and the ion migration barrier. Fig. 10 represents
memory endurance and retention temperature as a function of
calculated Emigration. The number of cycles reachable by a cell
decreases with Emigration [Fig. 10(a)], while temperature stabil-
ity increases [Fig. 10(b)]. This phenomenon can be understood
with a simple description of ion migration in the resistive
layer. Fig. 10 shows how endurance is facilitated by a low
migration barrier; a low energy will lead to an easy formation
and disruption of the ﬁlament. A high migration energy barrier
will lead to stable retention, making the ﬁlament harder to be
dissolved. In HfO
2,s i n c eC ui VO and VO calculated migration
energies are of the same order of magnitude (see Table III),
we considered a hybrid ﬁlament resulting from both Cu and
VO diffusion. Thus, Emigration is calculated as the average
value between Cu i VO and VO in Fig. 10. It clearly appears
that the ion migration energy is linked to the tradeoff between
memory endurance and retention.
Based on this observation, we propose the following
description: an RRAM device can be characterized by a
maximum energy E
max that can be sustained before oxide
degradation and breakdown. This maximum energy is reached
when a maximum number of cycles, Ncmax, is achieved.
The RRAM is then led to irreversible breakdown. For each
cycle, the provided energy Ecycle depends on the programming
conditions. The maximum energy can thus be expressed with
Ncmax and Ecycle [eq. 3(a)]. A fraction of Ecycle is transferred
to ions and serves to switch the RRAM ( Eswitch) while another
is transferred to the network and leads to material degradation
such as defect generation. In our case, the endurance failure
is detected when RRAM cannot be erased anymore. The cell
is then stuck in R
ON state, meaning that there are too much
defects in the oxide, which leads to a permanent conductive
path. Creation of defects occurs during cycling and can be
accelerated by the increase of V
RESET [18].
As it was previously demonstrated in Section III, increasing
ROFF /RON ratio leads to a drop in the number of cycles a
cell can sustain. Thus, the maximum energy is reached faster
when the WM increases. In our description, we link the WM
to the number of atoms involved during ﬁlament formation
and disruption ( Nmoved atoms ). Increasing SET programming
condition leads to bigger conductive path, composed of more
atoms and resulting in lower RON and higher WM. However,
increasing ROFF means a bigger free defect area in a resistive
layer. This means that the increasing RESET condition (time
and voltage) is equivalent to increasing the number of moving
atoms during RESET. On the other hand, we also demon-
strated that ion migration barrier impacts memory endurance.
If E
migration increases, the number of achieved cycle decreases.
In our description, we thus assume that Emigration is the
required energy to move and drift each ion in the oxide, during
ﬁlament formation and disruption. Switching energy can then
be expressed as E
switch = Emigration ∗ Nmoved atoms .T h e
maximum energy of the dielectric can thus be expressed with
the maximum number of cycles, the WM (through N moved
atoms) and E
migration. In [18], a similar description of this
phenomenon is given. Degradation of the oxide is controlled
by activation energy and V
RESET. Irreversible endurance failure
is associated to a high defect injection in conductive ﬁlament.
As explained previously, the switching energy (transferred
to ﬁlament ions) only represents a fraction of the total supplied
energy. Energy surplus is transferred to the network and
degrades the oxide. We can express the total energy provided
with Eswitch and a degradation parameter α: Eswitch = α ∗
Ecycle,0 <α< 1. α evolves with oxide degradation: the more
the RRAM is cycled, the more the oxide is degraded, leading
to a higher amount of energy transferred to the network.
We thus assume α = 1 for a fresh cell (without defect) and
gradually decreases with cycling. Finally we can express the
maximum energy that the dielectric can sustain depending
on N
cmax, Eswitch,a n d α [see eq. 3(b)] or depending on
endurance, retention, and WM [eq. (4)]
Emax ∝ Ncmax ∗ Ecycle (3a)
Ecycle ∝ Eswitch
α ∝ Emigration ∗ Nmoved atoms
α (3b)
Emax ∝ Ncmax ∗ Emigration ∗ Nmoved atoms
⇓⇓ ⇓
Endurance Retention WM (4)

---

<!-- PDF page 6 -->

This article has been accepted for inclusion in a future issue of this journal. Content is final as presented, with the exception of pagination.
6 IEEE TRANSACTIONS ON ELECTRON DEVICES
Fig. 11. (a) Migration energy barrier experimental and simulated
as a function of maximum endurance and (b) extraction example of
Eexperimental for W/HfO2).
where Emax the maximum energy a dielectric can sustain
before breakdown, Ncmax the maximum number of achieved
cycles, Ecycle the total energy supplied per cycle, Eswitch the
energy needed to move the a toms to create the ﬁlament,
α the degradation parameter, such as Eswitch = α∗ Ecycle and
0 <α< 1, Emigration the activation energy to move one atom
in the oxide and Nmoved atoms the number of involved atoms
in the ﬁlament formation and disruption.
In parallel, Emigration can be associated to the activa-
tion energy experimentally extracted from retention measure-
ments. Fig. 11(a) shows correlations between simulated and
experimental Emigration and memory endurance. Very good
agreement is observed between simulated and experimental
Emigration. Temperature stability failure can then be well asso-
ciated with the movement of Cu and VO atoms, which will
disrupt or form again the ﬁlament [27]. The ROFF failure,
in Fig. 11(b), reveals a ﬁlament formation.
In conclusion, the maximum dielectric energy is linked
to the memory performances and microscopic features of
the memory technology: WM (through the number of atoms
involved in the ﬁlament formation and disruption), mem-
ory endurance (maximum number of cycles), and memory
retention (through ion migration barrier). Consequently the
following hold.
1) Stronger programming conditions (current,
time, or voltage) enlarge the WM, decreasing N
cmax
(and endurance) for a given memory technology
(ﬁxed Emax), as observed in Fig. 6(a).
2) Changing memory stack, it makes possible to tune the
migration energy and favor either endurance or retention
[Figs. 6(b)and 10], for a ﬁxed WM.
The following observations must also be made.
1) Dielectric optimization is o f great importance; improv-
ing the material quality can lead to E
max increase and
the general improvement of the WM/endurance/retention
tradeoff. In the literature, J. Woo et al. [28] demon-
strated that the breakdown strength of a material can be
improved by thermal treatment.
2) TE/oxide combination needs to be optimized to tune
migration barrier. Unfortunately, this paper needs to be
done for each combination. However, some guidance
can be given to avoid high diffusion in oxides such as:
controlling oxide grain boundari es, material structure,
Fig. 12. Radar plot of the performance tradeoff for the various RRAM
classes reported in this paper, and critical key parameters (material
microscopic properties and operating conditions).
and material density. It can be added that small metal
atoms are an advantage for interstitial diffusion, which
is usually the most favorable.
3) Understanding defect generation during memory
endurance is also critical. For a given WM, different
sets of programming parameters (current, time, and
voltage) and memory stack features (TE and BE) should
minimize endurance failure, preventing early oxide
degradation ( Fig. 4 ). Moreover, smart programming
procedures are a way to adjust the energy provided to
the system, reducing the energy surplus transferred to
the network and degrading the memory dielectric.
Finally, the link between memory technology microscopic
features and memory performances can be used to design
an optimized memory stack. However, operating conditions
need to be adjusted to reach the speciﬁcations of a given
memory application. Fig. 12 summarizes optimization for
three applications. Low migration barrier and optimized WM
favor endurance for storage class memory (proposed stack:
W/HfO
2/CuTex ). Optimized SET current and long time pro-
gramming conditions enlarge the WM for data storage (pro-
posed stack: TiN/Al 2O3/CuTex ). High migration barrier and
large WM are favorable for embedded applications requiring
high stability (proposed stack: TiN/GdO
x /CuTex ).
VII. C ONCLUSION
In this paper, we demonstrated the link between microscopic
ﬁlament properties and memory WM, endurance, and retention
in the RRAM. In particular, we show that the ion migration
barrier controls the endurance/retention compromise. It can
be tuned by TE/oxide combinaison and material engineer-
ing. We quantiﬁed how the WM improvement degrades
endurance. The optimization of the dielectric quality is a way
to improve the overall WM/retention/endurance compromise.
Based on the understanding of the ﬁlament composition and
characteristics, various memory technologies are proposed to
reach the speciﬁcations of diff erent memory applications.
R
EFERENCES
[1] S. Sills et al. , “A copper ReRAM cell for storage class memory
applications,” in Proc. Symp. VLSI Technol. , Jun. 2014, pp. 1–2,
doi: 10.1109/VLSIT.2014.6894368.

---

<!-- PDF page 7 -->

This article has been accepted for inclusion in a future issue of this journal. Content is final as presented, with the exception of pagination.
NAIL et al.: LINK BETWEEN CBRAM PERFORMANCES AND MATERIAL MICROSCOPIC PROPERTIES 7
[2] J. R. Jameson et al. , “Conductive-bridge memory (CBRAM) with
excellent high-temperature retention,” in IEDM Tech. Dig. , Dec. 2013,
pp. 30.1.1–30.1.4, doi: 10.1109/IEDM.2013.6724721.
[3] Z. Wei et al., “Distribution projecting the reliability for 40 nm ReRAM
and beyond based on stochastic differential equation,” in IEDM Tech.
Dig., Dec. 2015, pp. 7.7.1–7.7.4, doi: 10.1109/IEDM.2015.7409650.
[4] ITRS. (2015). Executive Report . [Online]. Available: http://www.itrs2.
net/
[5] R. Strenz, “Embedded Flash tec hnologies and their applications:
Status & outlook,” in IEDM Tech. Dig. , Dec. 2011, pp. 9.4.1–9.4.4,
doi: 10.1109/IEDM.2011.6131521.
[6] C. Y . Chen, L. Goux, A. Fantini, A. Redolﬁ, G. Groeseneken, and
M. Jurczak, “Doped Gd-O based RRAM for embedded application,”
in Proc. IEEE 8th Int. Memory Workshop (IMW) , May 2016, pp. 1–4,
doi: 10.1109/IMW.2016.7495266.
[7] M. Barci et al. , “Impact of SET and RESET conditions on CBRAM
high temperature data retention,” in Proc. IEEE Int. Rel. Phys. Symp.,
Jun. 2014, pp. 5E.3.1–5E.3.4, doi: 10.1109/IRPS.2014.6860677.
[8] A. Fantini et al. , “Intrinsic switching variability in HfO 2 RRAM,”
in Proc. 5th IEEE Int. Memory Workshop , May 2013, pp. 30–33,
doi: 10.1109/IMW.2013.6582090.
[9] U. Celano et al. , “Three-dimensional observation of the conductive
ﬁlament in nanoscaled resistive memory devices,” Nano Lett. , vol. 14,
no. 5, pp. 2401–2406, May 2014, doi: 10.1021/nl500049g.
[10] S. Menzel, U. Böttger, M. Wimmer, and M. Salinga, “Physics of the
switching kinetics in resistive memories,” Adv. Funct. Mater. , vol. 25,
no. 40, pp. 6306–6325, Oct. 2015, doi: 10.1002/adfm.201500825.
[11] D. S. Jeong et al. , “Emerging memories: Resistive switching mecha-
nisms and current status,” Rep. Prog. Phys., vol. 75, no. 7, p. 076502,
2012, doi: 10.1088/0034-4885/75/7/076502.
[12] R. Waser, R. Dittmann, G. Staikov, and K. Szot, “Redox-based resis-
tive switching memories nanoionic mechanisms, prospects, and chal-
lenges,” Adv. Mater., vol. 21, nos. 25–26, pp. 2632–2663, Jul. 2009,
doi: 10.1002/adma.200900375.
[13] M. Saadi et al., “On the mechanisms of cation injection in conducting
bridge memories: The case of HfO 2 in contact with noble metal anodes
(Au, Cu, Ag),” J. Appl. Phys. , vol. 119, no. 11, p. 114501, Mar. 2016,
doi: 10.1063/1.4943776.
[14] G. Molas et al. , “Controlling oxygen vacancies in doped oxide based
CBRAM for improved memory performances,” in IEDM Tech. Dig. ,
Dec. 2014, pp. 6.1.1–6.1.4, doi: 10.1109/IEDM.2014.7046993.
[15] C. Nail et al. , “Understanding RRAM endurance, retention and win-
dow margin trade-off using experimental results and simulations,” in
IEDM Tech. Dig., Dec. 2016, pp. 4.5.1–4.5.4, doi: 10.1109/IEDM.2016.
7838346.
[16] C. Nail et al., “Atomistic mechanisms of copper ﬁlament formation and
composition in Al2o3-based conductive bridge random access memory,”
J. Appl. Phys., vol. 122, no. 2, p. 024503, Jul. 2017. [Online]. Available:
http://aip.scitation.org/doi/full/10.1063/1.4990979
[17] M. Barci et al. , “Bilayer metal-oxide conduc tive bridge memory tech-
nology for improved window margin and reliability,” IEEE J. Electron
Devices Soc., vol. 4, no. 5, pp. 314–320, Sep. 2016, doi: 10.1109/JEDS.
2016.2567219.
[18] S. Balatti et al. , “Pulsed cycling operation and endurance failure
of metal-oxide resistive (RRAM),” in IEDM Tech. Dig. , Dec. 2014,
pp. 14.3.1–14.3.4, doi: 10.1109/IEDM.2014.7047050.
[19] J. M. Soler et al., “The SIESTA method for ab initio order-N materials
simulation,” J. Phys., Condens. Matter , vol. 14, no. 11, p. 2745, 2002,
doi: 10.1088/0953-8984/14/11/302.
[20] P. Hohenberg and W. Kohn, “Inhomogeneous electron gas,” Phys. Rev.,
vol. 136, no. 3B, pp. B864–B871, Nov. 1964, doi: 10.1103/PhysRev.136.
B864.
[21] W. Kohn and L. J. Sham, “Self-cons istent equations including exchange
and correlation effects,” Phys. Rev., vol. 140, no. 4A, pp. A1133–A1138,
Nov. 1965, doi: 10.1103/PhysRev.140.A1133.
[22] N. Troullier and J. L. Martins, “Efﬁcient pseudopotentials for plane-wave
calculations,” Phys. Rev. B, Condens. Matter , vol. 43, pp. 1993–2006,
Jan. 1991, doi: 10.1103/PhysRevB.43.1993.
[23] L. G. Ferreira, M. Marques, and L. K. Teles, “Approximation to density
functional theory for the calculation of band gaps of semiconductors,”
Phys. Rev. B, Condens. Matter , vol. 78, no. 12, p. 125116, Sep. 2008,
doi: 10.1103/PhysRevB.78.125116.
[24] G. E. Delgado, A. J. Mora, M. Pirela, A. Velásquez-Velásquez,
M. Villarreal, and B. J. Fernández, “Structural reﬁnement of the
ternary chalcogenide compound Cu
2gete3 by X-ray powder diffraction,”
Phys. Status Solidi A , vol. 201, no. 13, pp. 2900–2904, Oct. 2004,
doi: 10.1002/pssa.200406850.
[25] W. Zheng, K. H. Bowen, Jr., J. Li, I. D¸ abkowska, and M. Gutowski,
“Electronic structure differences in ZrO 2 vs HfO2,” J. Phys. Chem. A ,
vol. 109, no. 50, pp. 11521–11525, Dec. 2005, doi: 10.1021/jp053593e.
[26] C. Freysoldt et al. , “First-principles calculations for point defects in
solids,” Rev. Mod. Phys., vol. 86, no. 1, pp. 253–305, Mar. 2014,
doi: 10.1103/RevModPhys.86.253.
[27] J. Guy et al. , “Guidance to reliability improvement in CBRAM using
advanced KMC modelling,” in Proc. IEEE Int. Rel. Phys. Symp. (IRPS) ,
Apr. 2017, pp. PM-2.1–PM-2.5, doi: 10.1109/IRPS.2017.7936384.
[28] J. Woo, A. Belmonte, A. Redolﬁ, H. Hwang, M. Jurczak, and L. Goux,
“Introduction of WO 3 layer in a Cu-based Al 2o3 conductive bridge
RAM system for robust cycling and large memory window,” IEEE
J. Electron Devices Soc. , vol. 4, no. 3, pp. 163–166, May 2016,
doi: 10.1109/JEDS.2016.2526632.
Authors’ photographs and biographies not available at the time of
publication.