### **Layer specific dendritic spine dynamics in primary motor cortex (M1) during learning**
**Nikhita Kaushik |  SURFiN Fellow @ Komiyama Lab, 8/25-5/26 | Mentor: Jennifer Li** 

---

#### **Project Overview**
Dendritic spines undergo structural changes during motor learning, but whether these changes differ across cortical layers remains unclear. Here, we used longitudinal two-photon microscopy to track individual spines on L2/3 and L5 dendrites in primary motor cortex (M1) across 14 consecutive days of training & imaging, comparing structural spine dynamics between layers. 

#### **Dataset**
- **Mice:** Wild-type mice (P60–P90), n=3
- **Cortical layers:** L2/3 & L5
- **Imaging duration:** 14 days, with water-restricted motor training
- **Imaging modality:** Two-photon microscopy at 1000 nm; red + green channels acquired simultaneously
- **Volumes acquired:** 1 µm steps, 20–60 µm depth
- **Viral labeling:** L2/3: CamKII-Cre + FLEX-tdTomato (200–300 µm depth); L5 corticospinal neurons: retrograde Cre + FLEX-eGFP (injected C4–C6)
- **Spine annotation:** Custom GUI; volume measurements validated with correlated electron microscopy

#### **Experimental Timeline**
- **Week 0** — Cranial window surgery
- **Week 2** — Begin water deprivation
- **Week 4–6** — Daily two-photon imaging + lever-press training (14 days)

#### **Analysis Pipeline**
Spines are tracked across all 14 imaging days per FOV, the following is measured:
- **Spine density** — spines per µm over time, raw & day-1 normalised
- **Spine turnover** — daily addition & elimination rates per µm
- **Spine survival** — fraction of pre-existing & newly formed spines surviving to later days
- **Spine lifetime** — distribution of how long individual spines persist
- **Spine volume change** — pairwise volume changes classified by plasticity type
- **Plasticity counts** — number of spines per plasticity category per dendrite per day
- **Volume × lifetime** — relationship between initial spine size & spine lifetime duration

All quantifications are compared between L2/3 & L5 using a mixed effects model (layer x day interaction, mouse as random effect) to control for repeated measurements within animals. 
