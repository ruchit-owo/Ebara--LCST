# Sources

SMILES rule: each comonomer is written as the neutral free monomer (vinyl form, acids protonated). For DMMDO_deprotected the repeat unit only exists after post-polymerisation acetal removal, so its SMILES is the enol monomer that would give the same repeat unit.

Search route: Europe PMC open-access full text (queries on NIPAM/NIPAAm/isopropylacrylamide + LCST/cloud point + copolymer), tables with LCST/Tcp columns pulled from the XML and checked by hand against the methods section.

## Training set (role = train)

| source key | citation | DOI | used |
|---|---|---|---|
| Chen2025_Materials_18_1614 | Chen Z et al., Materials 2025, 18, 1614 | 10.3390/ma18071614 | Table 2, DHPMA |
| Bucatariu2024_Gels_10_834 | Bucatariu S et al., Gels 2024, 10, 834 | 10.3390/gels10120834 | Table 1, HEAAm |
| Abdelaty2023_ACSOmega_8_45026 | Abdelaty MSA et al., ACS Omega 2023, 8, 45026 | 10.1021/acsomega.3c06634 | Table 1 groups A-C (HEMA, DMAAm, styrene); terpolymers not used |
| Kashima2025_Macromolecules_58_8753 | Kashima R et al., Macromolecules 2025, 58, 8753 | 10.1021/acs.macromol.5c01401 | cloud point table, DMMDO before/after deprotection |
| Safakas2021_IJMS_22_3824 | Safakas K et al., Int J Mol Sci 2021, 22, 3824 | 10.3390/ijms22083824 | Table 1, NtBAm (free grafting chains) |
| Chen2018_Polymers_10_1259 | Chen BY et al., Polymers 2018, 10, 1259 | 10.3390/polym10111259 | Table 1, HMAAm (NMA) with ~0.4% rhodamine comonomer |
| Liang2017_Polymers_9_136 | Liang FC et al., Polymers 2017, 9, 136 | 10.3390/polym9040136 | Table 1 P1 only (P2 is a terpolymer) |
| Fundueanu2022_Pharmaceutics_14_865 | Fundueanu G et al., Pharmaceutics 2022, 14, 865 | 10.3390/pharmaceutics14040865 | Table 1, VI, water column |
| Li2019_RSCAdv_10_323 | Li J et al., RSC Adv 2019, 10, 323 | 10.1039/c9ra08832f | Table 1, AAc at pH 5.7 and 6.6 |
| Nagase2019_IJMS_20_430 | Nagase K et al., Int J Mol Sci 2019, 20, 430 | 10.3390/ijms20020430 | Table 1 + text, DMAPAAm and DMAAm |
| Nieswandt2021_Polymers_14_62 | Nieswandt K et al., Polymers 2021, 14, 62 | 10.3390/polym14010062 | Table 3, DMAAm |
| Tugaeva2026_ACSOmega_11_22663 | Tugaeva GK et al., ACS Omega 2026, 11, 22663 | 10.1021/acsomega.5c11071 | Table 1, NtBA statistical copolymer only |

## Ebara group check set (role = ebara_check, not for training)

| source key | citation | DOI | used |
|---|---|---|---|
| Techawanitchai2011_STAM_12_044609 | Techawanitchai P, Yamamoto K, Ebara M, Aoyagi T, Sci Technol Adv Mater 2011, 12, 044609 | 10.1088/1468-6996/12/4/044609 | N90H10 / N80H20 / N70H30 values; read from PMC page summary, VERIFY against PDF |
| Niiyama2018_Polymers_10_1018 | Niiyama E, Uto K, Lee CM, Sakura K, Ebara M, Polymers 2018, 10, 1018 | 10.3390/polym10091018 | 21 mol% HMAAm, 43 C |
| Chen2022_FrontBioeng_10_1046147 | Chen L, Fujisawa N, Takanohashi M, Ebara M, Front Bioeng Biotechnol 2022, 10, 1046147 | 10.3389/fbioe.2022.1046147 | 20 mol% HMAAm, 42.1 C |

Leakage flag: HMAAm (NMA) also appears in the training set (Chen2018, Liang2017). For the Ebara check those rows have to be dropped or the check run as a held-out-comonomer fold.

## Looked at, not used

| paper | reason |
|---|---|
| Ebara group, Polymers 2019, 11, 42 (PMC6401969), NIPAAm-co-PBA | only a 7-14 C range from an SI figure |
| Ebara group, STAM 2012, 13, 064203 (dextran nanofibers) | publisher blocked, not in PMC full text |
| PMC10658503, NIPAM-co-SBMA | polymerised in situ inside a gel; LCST method not stated |
| PMC6630648, NIPAM-co-dopamine methacrylamide | values only in figures, 4-15 wt% solutions |
| PMC8473767, NIPAM-HEMA by ultrasound | made by chain scission + grafting, not a random copolymer |
| PMC8382249, NIPAM-co-mannosyl methacrylates | turbidity onset flat, LCST only by NMR |
| PMC9143634, PMC12783973, PMC11679721 | crosslinked microgels (VPTT) |
| PMC5349575 | copolymers grafted onto graphene oxide |
| PMC4025585, PMC3208240, PMC3446841, PMC7023078, PMC9611013 | terpolymers |
| PMC10048633 | same NtBAM chains as Safakas2021 (duplicate) |
| PMC8912701 (review table, AAm 10-30%) | secondary source; original not open |
