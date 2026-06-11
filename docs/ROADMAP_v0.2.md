# SpatialHAC.jl — Roadmap d'implémentation v0.2 *(vérifiée 2026-06-10)*

> **But** : faire passer le package de prototype mono-fonction (v0.1.0, `vcov_conley` + `suggest_cutoff`, 368 lignes) à une bibliothèque autonome, *feature-complete*, publiable dans *Journal of Open Source Software*.
>
> **Principe non négociable** : *chaque pièce est validée deux fois (définition dense interne + référence externe indépendante), property-testée, auditée de façon adversariale avant d'être taggée.* Toute fonctionnalité non validable contre une référence externe est **coupée** (JOSS : « no half-baked solutions »).
>
> ⚠️ **Ce document a été vérifié** : toutes les références externes et les claims structurels ont été contrôlés contre les sources le 2026-06-10. Voir le **§ Journal de vérification** en fin de fichier pour le verdict claim-par-claim et les corrections appliquées par rapport au premier jet.

---

## 0. Pré-requis bloquants découverts à la vérification *(à corriger AVANT tout tag)*

Ces points ne sont pas des « nouvelles features » : ce sont des défauts du package **actuel** révélés par la vérification. Ils doivent être réglés avant la v0.2.0.

### 0.1 — Export fantôme `vcov_cluster` / `ClusterResult` *(CONFIRMÉ localement)*
La ligne 36-37 de `src/SpatialHAC.jl` **exporte** `vcov_cluster` et `ClusterResult`, mais `isdefined(SpatialHAC, :vcov_cluster) == false` (vérifié par chargement du module). La suite passe uniquement parce qu'aucun test n'y touche. **`Aqua.test_all` lèvera un échec d'export indéfini** (chantier 9). → Soit on implémente le chantier 2 (qui définit ces symboles), soit on retire l'export. Comme le chantier 2 est retenu, l'export devient valide une fois l'implémentation faite ; entre-temps c'est une dette à ne pas oublier.

### 0.2 — Claim « Bartlett PSD-admissible » SURÉVALUÉ dans le package livré *(CONFIRMÉ par lecture Kelejian-Prucha 2007)*
Le README (l. 29-31), le docstring `src` et `docs/src/estimator.md` affirment :
> *Bartlett kernel `K(d)=max(0,1−d/cutoff)` (PSD-admissible; Kelejian & Prucha 2007)*

**C'est inexact en 2-D.** Kelejian-Prucha (2007, Remark 1, via Golubov 1981) montrent que la famille triangulaire `Kₙ(x)=(1−x)ⁿ₊` appartient à la classe de Schoenberg `Pₚ` **ssi `n ≥ (p+1)/2`**. Le Bartlett linéaire `K₁` (= notre noyau) est donc dans `P₁` **seulement** : PSD garanti en **1-D**, pas en 2-D géographique (qui exige `n≥2`, le noyau triangulaire au carré `K₂`). KP2007 avertissent eux-mêmes que l'estimateur « may not be positive semi-definite in finite samples ».
**Conséquence pour le package existant** : il faut (a) **corriger le wording** README/src/docs — ne plus écrire « PSD-admissible » sans qualificatif, et présenter le **flooring de valeurs propres** comme le vrai filet de sécurité PSD (déjà implémenté) ; (b) le chantier 3 doit ajouter `K₂` (triangulaire au carré) comme le noyau réellement PSD-garanti en 2-D. Reformulation correcte recommandée :
> *Bartlett `K₁` est PSD-garanti en 1-D (Schoenberg class `P₁`) ; en 2-D, aucun noyau triangulaire linéaire n'est PSD-garanti — utiliser `K₂` (∈`P₂`) ou s'appuyer sur le flooring de valeurs propres. Conley (1999) et KP2007 ancrent la PSD dans la positive-définitude de Bochner/Schoenberg.*

### 0.3 — Manifest périmé *(CONFIRMÉ : warning Pkg au chargement)*
`Project.toml` a gagné la dép `Random` sans `Pkg.resolve()` → warning « project dependencies or compat requirements have changed since the manifest was last resolved ». Lancer `Pkg.resolve()` et committer le `Manifest.toml` à jour avant le tag.

---

## Séquencement de release (contrainte ferme)

```
v0.1.0 (en cours d'enregistrement General Registry)
   │  PR JuliaRegistrator auto-merge ~2026-06-13 → TagBot tag v0.1.0
   ▼
v0.2.0  ← regrouper TOUS les chantiers retenus + pré-requis §0, PUIS :
   │  1. Pkg.resolve() ; bump Project.toml version → 0.2.0
   │  2. @JuliaRegistrator register (commentaire sur le commit de release)
   │  3. attendre tag v0.2.0 ; archiver (Zenodo → DOI)
   ▼
Soumission JOSS sur le tag v0.2.0 (paper.md + paper.bib)
```

**Note JOSS vérifiée** : il n'existe **plus** de seuil dur « 1000 lignes » dans les critères JOSS actuels (révision 2026). Le vrai barème est qualitatif : *substantial scholarly effort*, ~3 mois de travail individuel minimum, *feature-complete (no half-baked solutions)*, historique public idéalement ≥6 mois, « likely to be cited ». → **Stratégie : ne pas gonfler le nombre de lignes ; soigner la complétude, la nouveauté (niche vide, vérifiée), les tests et l'historique.** La cible « ~1000 lignes » reste un repère de maturité, pas un gate.

---

## Vue d'ensemble des chantiers

| # | Chantier | Validation externe | Priorité | Bloquant JOSS |
|---|----------|--------------------|----------|---------------|
| 0 | **Pré-requis** : export fantôme, wording PSD, Manifest | — (corrections internes) | **P0** | Oui |
| 1 | Support modèles pondérés | sandwich dense pondéré + R `clubSandwich` (lmer pondéré) | **Haute** | Oui (limitation affichée) |
| 2 | `vcov_cluster` (CR0 / CR1 / CR1S) | R `clubSandwich::vcovCR` types `"CR0"`,`"CR1"`,`"CR1S"` | **Haute** | Oui (2ᵉ fonctionnalité) |
| 3 | Noyaux (`K₂` triangulaire², uniforme, Epanechnikov) + cadre Schoenberg | covariance dense + test PSD Monte Carlo | Moyenne | Non |
| 3.5 | **Diagnostics spatiaux exposés** : `covariogram` Ĉ(h) + `variogram`/`semivariogram` γ̂(h) | brute-force dense + R `gstat::variogram` + `GeoStats.jl` | Moyenne | Oui (3ᵉ capacité visible) |
| 4 | Distances projetées (euclidien/planaire) | distance euclidienne dense vs haversine | Moyenne | Non |
| 5 | API `StatsAPI` / `show` / `CoefTable` | comparaison vs `coeftable(m)` | Moyenne | Oui (ergonomie) |
| 6 | HAC spatio-temporel | **SEULEMENT si** validable vs Driscoll-Kraay (`plm::vcovSCC`) ; sinon COUPER | Basse | Non |
| 7 | Vignette de validation (docs exécutables) | reproduit les tableaux de validation | Moyenne | Oui (article) |
| 8 | `paper.md` (750-1750 mots) + `paper.bib` (BibLaTeX) | — | — | Oui |
| 9 | Qualité : `Aqua.jl`, `JET.jl`, docstrings, exemples reproductibles | `Aqua.test_all`, `JET.test_package` | Moyenne | Oui |

---

## Chantier 1 — Support des modèles pondérés *(ACTIF, math confirmée)*

### Pourquoi
Modèles mixtes pondérés très fréquents : pondération par surface (area-weighted), effort d'échantillonnage, inverse-variance, correction de biais. Actuellement `vcov_conley` lève une `ArgumentError` dès qu'il y a des poids (l. 120-122).

### Mathématiques *(confirmée 2 façons : probe empirique + lecture internals)*
Dans l'espace blanchi par les poids, mettre X, ê, W à l'échelle par `√w` ; machinerie identique ensuite.
```
sqw = isempty(m.sqrtwts) ? ones(n) : Float64.(m.sqrtwts)   # m.sqrtwts = √(prior weights)
X̃ = sqw .* X ;  ẽ = sqw .* ê ;  W̃ = sqw .* W              # ê = response(m) − X·β̂ ; W = ZΛ brut
```
**Vérifié localement (2026-06-10)** : sur un LMM pondéré, `m.sqrtwts` non-vide et `≈ sqrt.(w)` ; vide pour un fit non pondéré ; `response(m)` = y original ; `scaled_re_matrix` renvoie ZΛ **brut** (non pondéré) → doit être scalé. **Probe (session antérieure)** : self-check `varest(m)·B ≈ vcov(m)` → 7.4e-14 avec ce scaling (Interp. A), 0.78 / 1.6e-3 pour les variantes fausses. Le garde-fou runtime existant attrape donc tout double-scaling.

### Travail
1. `src` : remplacer le bloc garde (l. 120-122) par `sqw = isempty(m.sqrtwts) ? ones(n) : Float64.(m.sqrtwts)` ; appliquer `Xw/ehat_w/W_w` dans F/OinvX/bread/S. (`n` déjà défini via `size(X,1)`.)
2. `suggest_cutoff` : conserver ê **marginal brut** (Lehner travaille sur résidus non blanchis) ; ajouter une phrase docstring pour les modèles pondérés. Pas de changement de code.
3. Test : « weighted == dense weighted sandwich » — fitter `wts=df.w`, étendre `dense_sandwich` (scaler les lignes par √w), concordance < 1e-12.
4. **Validation externe** : `test/crosscheck_weighted.R` — `lmer` pondéré + sandwich Conley pondéré manuel, concordance ~1 % à θ apparié. (clubSandwich gère `lmerMod` — confirmé — utilisable comme 2ᵉ ancre cluster, mais pour le Conley *spatial* l'ancre reste le sandwich manuel R, comme `crosscheck_conley_lme4.R`.)
5. Docs : retirer « Unweighted models only » du caveat 3 README ; exemple pondéré dans `estimator.md`.

**Fini** = self-check runtime OK pondéré + dense < 1e-12 + crosscheck R < 1 % + README à jour + suite verte.

---

## Chantier 2 — `vcov_cluster` : SE cluster-robustes (CR0 / CR1 / CR1S)

> ⚠️ **Correction vérifiée** : la nomenclature CR doit être exacte. D'après le manuel `clubSandwich` v0.7.0, les types sont `CR0, CR1, CR1p, CR1S, CR2, CR3`. **CR1 ≠ ce que j'avais écrit** : CR1 = CR0 × `m/(m−1)` seulement. Le facteur complet `(m(N−1))/[(m−1)(N−p)]` que je nommais « CR1 » est en réalité **CR1S** (le défaut Stata). On implémente et on nomme correctement : **CR0**, **CR1** (`m/(m−1)`), **CR1S** (facteur complet). CR2/CR3 hors scope (corrections type-HC2/HC3 par cluster, infaisables à 7,8 M lignes, inutiles avec des milliers de clusters).

### Mathématiques
Même bread GLS `B=(X'Ω⁻¹X)⁻¹` ; viande clusterisée par sommation de blocs (jamais d'inversion globale) :
```
M_CR0 = Σ_g (Σ_{i∈g} rᵢ)(Σ_{i∈g} rᵢ)'      rᵢ = (Ω⁻¹X)ᵢ·êᵢ
CR1  = (m/(m−1))·M_CR0
CR1S = (m(N−1))/[(m−1)(N−p)]·M_CR0          (= défaut Stata)
V = B · M · B
```
CR0 = sandwich Liang-Zeger 1986 (« the original form ... no small-sample correction » — manuel clubSandwich). Port du code moteur validé `scripts/julia/fire_mixed_models/conley_se.jl`.

### Travail
1. `vcov_cluster(m, cluster_id; type=:CR1)` → `ClusterResult` (`vcov, se, n_clusters, type, dof`). **Ceci résout aussi le §0.1** (l'export devient défini).
2. Multi-way (Cameron-Gelbach-Miller 2011, `V = V₁ + V₂ − V₁₂` — **formule confirmée**) *si* validable ; sinon one-way seulement.
3. Réutiliser le garde-fou runtime.
4. Tests : (i) dense par blocs interne < 1e-12 ; (ii) **cross-langage R `clubSandwich::vcovCR(model, type=...)`** pour `"CR0"`, `"CR1"`, `"CR1S"` à θ apparié < 1 % ; (iii) un cluster unique == HC0 ; (iv) limite OLS.
5. Docs + entrée tableau de validation, nommage CR exact.

**Fini** = concordance R clubSandwich CR0/CR1/CR1S < 1 % + dense < 1e-12 + scalabilité (pas d'allocation O(N²)).

---

## Chantier 3 — Noyaux additionnels + cadre Schoenberg PSD

> ⚠️ **Reframe vérifié** (voir §0.2). Le récit « Bartlett PSD vs uniforme non-PSD » est **faux en 2-D** : ni l'uniforme ni le Bartlett linéaire `K₁` ne sont PSD-garantis en 2-D. La PSD exige un noyau de la classe de Schoenberg `Pₚ`.

### Travail
1. `kernel::Symbol = :bartlett` (défaut, `K₁`) ; ajouter :
   - `:bartlett2` / `:triangular2` = `K₂(x)=(1−x)²₊` → **PSD-garanti en 2-D** (`∈P₂`, Golubov 1981) — le noyau honnête pour données spatiales.
   - `:uniform` (troncature dure, Conley 1999 original) — **non PSD-garanti**, fourni avec avertissement.
   - `:epanechnikov`.
2. **Garde PSD honnête** : conserver le flooring de valeurs propres comme filet pratique ; `@warn` explicite quand le noyau choisi n'est pas dans `Pₚ` pour la dimension utilisée ET `min_eig < −tol`. Documenter `K₂` (ou quadratic-spectral/gaussien) comme seuls choix PSD-garantis en 2-D ; ancrer dans Schoenberg/Bochner (Conley 1999 ; KP2007 fn. 16).
3. Tests : (i) chaque noyau vs covariance dense définitionnelle ; (ii) Monte Carlo PSD — fraction de réplicats `min_eig<−tol` par noyau (`K₂`≈0 en 2-D, `K₁`/uniforme >0 attendu → test de régression du warning) ; (iii) couverture comparée.
4. **Corriger en même temps le wording PSD du package livré** (§0.2) : README, src docstring, `estimator.md`.

**Fini** = chaque noyau == dense ; warning PSD testé sur cas dégénéré ; doc Schoenberg correcte ; wording livré corrigé.

---

## Chantier 3.5 — Diagnostics spatiaux exposés : covariogramme **et** variogramme

> **Pourquoi maintenant** : le covariogramme `Ĉ(h)` est **déjà calculé** mais enfoui dans `suggest_cutoff` (renvoyé dans `CovariogramResult.bins`/`.C`, jamais appelable seul). Le variogramme `γ̂(h)` — l'objet géostatistique que tout chercheur spatial reconnaît — **n'existe pas**. Les exposer transforme une mécanique interne en un **diagnostic de premier plan** : l'utilisateur veut *voir* la structure de corrélation spatiale de ses résidus (portée, sill, nugget) avant de choisir un cutoff, pas seulement recevoir un nombre. C'est la 3ᵉ capacité visible du package (à côté des deux estimateurs de SE) — bon pour le récit de complétude JOSS, et quasi-gratuit (même boucle de paires que le covariogramme déjà codé).

### Définitions (terminologie vérifiée — voir Journal)
- **Covariogramme** `Ĉ(h)` = moyenne par bin des produits `êᵢ·êⱼ` de paires à distance ≈ h (même période). C'est la fonction de covariance empirique des résidus. C'est *exactement* ce que `suggest_cutoff` calcule (méthode Lehner 2026).
- **Semivariogramme** `γ̂(h)` = `½ · moyenne[(êᵢ−êⱼ)²]` par bin. **`γ` est le SEMIvariogramme ; le « variogramme » au sens strict est `2γ`** (convention Cressie). Dans la pratique logicielle, ce que `gstat` et `GeoStats.jl` appellent « variogram » et renvoient est `γ` (le semivariogramme) — **donc pas de facteur 2 à appliquer** quand on les compare.
- **Identité de liaison** (sous stationnarité d'ordre 2, *sill fini*) : `γ(h) = Ĉ(0) − Ĉ(h)`, avec `Ĉ(0) = sill = variance`. Sert de **property-test interne croisant les deux fonctions**.

### Travail
1. **`covariogram(model_or_ehat, lat, lon, period; nbins=150, max_frac=2/3, max_points, rng, distance=:haversine)`** → expose le calcul interne déjà existant comme fonction publique. Renvoie un `SpatialDiagnostic` (champs `h` = centres de bins, `value` = Ĉ(h), `n_pairs`, `kind=:covariogram`).
2. **`variogram(model_or_ehat, lat, lon, period; estimator=:matheron, ...)`** → semivariogramme `γ̂(h)`. Estimateurs `:matheron` (½ des différences au carré) et `:cressie` (robuste, 4ᵉ puissance) — mêmes noms que `gstat`/`GeoStats.jl` pour la lisibilité. Même structure `SpatialDiagnostic` (`kind=:semivariogram`). Documenter clairement « renvoie γ (semivariogramme), pas 2γ ».
3. **Refactor DRY** : extraire la boucle de binning de paires en un helper interne `_binned_pairs(stat, ...)` paramétré par la statistique (`product` → covariogramme, `halfsqdiff` → semivariogramme). `suggest_cutoff` appelle alors `covariogram` en interne. Réduit la duplication et garantit que les trois (suggest_cutoff, covariogram, variogram) partagent exactement le même binning/sous-échantillonnage/convention même-période.
4. **Property-test de liaison** : `γ̂(h) ≈ Ĉ(0) − Ĉ(h)`. ⚠️ **Gotcha vérifié** : l'identité ne tient que si le sill est atteint (stationnarité d'ordre 2, `Ĉ(0)` fini). Le test doit **échouer bruyamment** si `max_frac·dmax` n'atteint pas un plateau (sinon désaccord silencieux). Reconstruire `Ĉ(0)` comme variance des résidus.
5. **Recette de tracé** : un exemple docs (ou recette Plots/Makie légère, dépendance optionnelle via package extension pour ne pas alourdir) traçant Ĉ(h) et γ̂(h) avec repères nugget/sill/range.
6. Exporter `covariogram`, `variogram`, `SpatialDiagnostic` ; déprécier proprement l'accès via `CovariogramResult` (ou le conserver, `suggest_cutoff` gardant son type de retour spécialisé qui embarque `cutoff`/`crossed`).

### Validation *(2 ancres externes confirmées)*
1. **Brute-force dense interne** < 1e-12 : double boucle naïve sur toutes les paires même-période, pour Ĉ(h) ET γ̂(h).
2. **Cross-package R `gstat`** : `variogram(z~1, data)` → colonne `gamma` = γ̂(h) ; `variogram(z~1, data, covariogram=TRUE)` → Ĉ(h) (⚠️ **gotcha vérifié** : la colonne reste nommée `gamma` même en mode covariogramme — s'assurer du flag, pas du nom de colonne). Aligner `cutoff`/`width` (gstat) ↔ `max_frac·dmax`/`nbins` (nous). Concordance < 1 %.
3. **Cross-package Julia `GeoStats.jl`** : `EmpiricalVariogram(data, var; nlags, maxlag, estimator=:matheron)`, accès `h, g, n = values(γ)`. ⚠️ **gotchas vérifiés** : `nbins` a été renommé `nlags` (v0.7.0) ; **GeoStats.jl n'a PAS de type covariance empirique** → il ne valide que le semivariogramme ; pour Ĉ(h) l'ancre externe reste `gstat` (`covariogram=TRUE`) ou la reconstruction `Ĉ(0)−γ̂` interne.
4. Property-test de liaison (point 4 ci-dessus) sur un champ GP simulé à sill fini.

**Fini** = `covariogram`/`variogram` == brute-force dense < 1e-12 ; == `gstat` (γ et C) < 1 % ; == `GeoStats.jl` (γ) < 1 % ; identité `γ=Ĉ(0)−Ĉ(h)` testée avec garde-plateau ; `suggest_cutoff` refactoré pour appeler `covariogram` (non-régression : mêmes résultats qu'avant) ; exemple de tracé dans les docs.

---

## Chantier 4 — Distances projetées (planaire / euclidien)

### Pourquoi
Haversine suppose lat/lon en degrés. Utilisateurs avec coords projetées (UTM/Albers, mètres) → haversine serait faux.

### Travail
1. `distance::Symbol = :haversine` (défaut) ; ajouter `:euclidean` (x/y planaires, même unité que `cutoff`).
2. Validation : euclidien dense vs grille ; sur petit domaine projeté localement, haversine≈euclidien (sanity).
3. Docs : `:euclidean` → `cutoff` dans l'unité des coords ; `:haversine` → degrés + `cutoff` km.

**Fini** = euclidien == dense ; doc unités sans ambiguïté ; test de cohérence.

---

## Chantier 5 — `StatsAPI` / `show` / `CoefTable`

> ⚠️ **Correction vérifiée** : `CoefTable` vit dans **StatsBase.jl** (pas StatsAPI), et son constructeur est **positionnel**, pas par champs nommés : `CoefTable(cols::Vector, colnms::Vector, rownms::Vector, pvalcol::Int=0, teststatcol::Int=0)`. `coeftable` (générique déclaré dans StatsAPI) doit renvoyer un `CoefTable` compatible Tables.jl.

### Travail
1. `StatsAPI.vcov`, `StatsAPI.stderror` (défaut `sqrt.(diag(vcov))`), `StatsAPI.coeftable` pour `ConleyResult`/`ClusterResult`.
2. `coeftable` construit positionnellement : `cols=[est, se, z, p, lo, hi]`, `colnms=["Coef.","Std. Error","z","Pr(>|z|)","Lower 95%","Upper 95%"]`, `rownms=coefnames`, `pvalcol=4`, `teststatcol=3`.
3. `Base.show(io, ::MIME"text/plain", r)` → tableau aligné (réutiliser le rendu CoefTable).
4. Stocker `coefnames` du modèle source dans le résultat.
5. Tests : estimations == `coef(m)` ; SE/z/p cohérents avec `r.se` ; `show` ne lève pas, contient les colonnes.

**Fini** = `coeftable(res)` produit un `CoefTable` valide (constructeur positionnel) ; `show` lisible ; estimations inchangées.

---

## Chantier 6 — HAC spatio-temporel — **COUPÉ (décision Thierry 2026-06-11)**

> **Statut : coupé, documenté comme limitation assumée** (README Caveats §1).
> Raison : pas d'ancre externe praticable — `plm::vcovSCC` (Driscoll-Kraay)
> cible les modèles panel `plm`, pas l'estimand GLS du modèle mixte ; livrer
> sans validation à deux voies indépendantes = exactement le risque du piège v1.
> Alternative documentée : `vcov_cluster` à niveau grossier pour sonder la
> dépendance inter-périodes. Le texte original du chantier est conservé
> ci-dessous pour référence si une ancre devient disponible.

### (archivé) Chantier 6 — HAC spatio-temporel *(conditionnel — couper si non validable)*

> Validation externe disponible et **confirmée** : `plm::vcovSCC` = estimateur Driscoll-Kraay (1998) « consistent with cross-sectional and serial correlation ». C'est l'ancre. Pas d'ancre → on coupe (JOSS « no half-baked »).

### Travail (si feu vert validation)
1. Noyau produit `K(d_space)·K_time(|tᵢ−tⱼ|)`, paires sur fenêtre `lag_max`.
2. Validation : (i) dense espace-temps interne ; (ii) `lag=0` == Conley actuel ; (iii) spatial-cutoff→∞ uniforme + temps == `plm::vcovSCC` (Driscoll-Kraay) à θ apparié < 1 %.
3. Si (iii) échoue après investigation → **retirer**, documenter pourquoi dans Caveats.

**Fini (ou coupe)** = concordance Driscoll-Kraay < 1 % → livrer ; sinon couper proprement.

---

## Chantier 7 — Vignette de validation (docs exécutables)

`docs/src/validation.md` avec blocs `@example`/doctest reproduisant : limite OLS == Conley brute-force ; récupération portée covariogram (GP sphérique) ; CR0 == clubSandwich ; (si livré) espace-temps == Driscoll-Kraay. Chaque ligne du tableau de validation README → un test reproductible. C'est l'argument central de l'article.

---

## Chantier 8 — Soumission JOSS

> **Exigences vérifiées (docs JOSS actuelles, rév. 2026)** :

`paper.md`, **750-1750 mots** (pas 250-1000), front-matter YAML : `title`, `tags`, `authors` (`name`, `orcid`, `affiliation`), `affiliations` (`name`, `index`), `date`, `bibliography: paper.bib`. Sections requises actuelles : **Summary** (audience non-spécialiste), **Statement of need**, **State of the field**, **Software design**, **Research impact**, **+ divulgation d'usage d'IA**. `paper.bib` en **BibLaTeX** (BibTeX accepté), citations `@auteur:année`.

**Gate dur d'acceptation** : licence OSI (✓ MIT), docs, tests, contrôle de version, **release taggée + archive DOI (Zenodo)**. Enregistrement General Registry = norme communautaire Julia (attendue) mais **pas** un gate JOSS strict ; le gate est le DOI archivé.

Positionnement (argument de nouveauté, niche vérifiée vide) : `conleyreg`/`fixest`/`acreg` = OLS/GLM ; `clubSandwich` = mixte mais **cluster-robuste only, non spatial** (confirmé) ; `plm::vcovSCC` = spatial-HAC mais sur modèles `plm`, **pas** mixtes. Le créneau « spatial-HAC Conley sur effets fixes d'un modèle mixte à l'échelle » est réellement vide.

`paper.bib` : Conley 1999 ; Liang-Zeger 1986 ; Kelejian-Prucha 2007 ; Golubov 1981 (PSD `Pₚ`) ; Cameron-Gelbach-Miller 2011 + Cameron-Miller 2015 ; Lehner 2026 ; **Huang, Wiedermann & Zhang 2022, MBR** (PAS « Huang 2022 » seul) ; Driscoll-Kraay 1998 (si ch. 6) ; Cressie 1993 (variogramme, convention γ/2γ — ch. 3.5) ; MixedModels.jl ; Julia.

---

## Chantier 9 — Qualité logicielle

1. `Aqua.test_all(SpatialHAC)` — **attrapera l'export fantôme §0.1 s'il reste** ; méthodes ambiguës, deps inutilisées, piracy.
2. `JET.test_package` — analyse statique des erreurs de type.
3. Docstrings complètes sur tout exporté (params, retour, exemple, refs).
4. Exemple reproductible README/`examples/` avec **données simulées** (pas de dép aux données m31n privées).
5. Couverture >90 % (badge Codecov).
6. `Pkg.resolve()` + Manifest committé (§0.3).

---

## Protocole de validation transversal *(s'applique à TOUTE pièce)*

1. **Définitionnelle dense** — formule sans Woodbury/grille, petit échantillon, ~1e-12.
2. **Externe indépendante** — autre langage/package (R surtout) à θ apparié, < 1 %.
3. **Property test** — invariant (limite OLS, cas dégénéré → estimateur connu, récupération de paramètre injecté).
4. **Audit adversarial** — re-dériver l'algèbre ; chercher le cas où la formule plausible diffère de la vraie (leçon v1 : Ω⁻¹-sur-résidu vs design coïncident à Ω=I → seuls les tests qui ancrent la FORMULE l'attrapent).
5. **Garde-fou runtime** — étendre `varest(m)·bread ≈ vcov(m)`.

> Échec de (2) sans explication → **coupé**. Un package plus petit et 100 % solide vaut mieux qu'un large et douteux.

---

## Ordre d'exécution recommandé

```
0. Pré-requis §0 (export fantôme via ch.2, wording PSD via ch.3, Pkg.resolve)
1. Chantier 1 (pondéré)      ← actif, math confirmée 2 façons
2. Chantier 2 (cluster CR0/CR1/CR1S) ← définit l'export fantôme, 2ᵉ feature
3. Chantier 5 (StatsAPI/show)← ergonomie (constructeur CoefTable positionnel)
4. Chantier 3 (noyaux + reframe PSD) ← corrige aussi le wording livré
4b. Chantier 3.5 (covariogram + variogram) ← refactor DRY de suggest_cutoff, 3ᵉ capacité visible
5. Chantier 4 (distances)    ← additif
6. Chantier 6 (espace-temps) ← CONDITIONNEL, valider vs plm::vcovSCC d'abord
7. Chantier 9 (Aqua/JET/docstrings/Manifest)
8. Chantier 7 (vignette) + 8 (paper.md 750-1750 mots)
9. Pkg.resolve → bump 0.2.0 → @JuliaRegistrator → tag → Zenodo DOI → JOSS
```

Chantiers 3/4/5 indépendants, parallélisables. Chantier 6 = seul à risque de coupe. Ne pas tagger v0.2.0 tant que (a) v0.1.0 enregistrée, (b) pré-requis §0 réglés, (c) tous les chantiers retenus passent le protocole, (d) suite + Aqua + JET vertes 3 OS.

---

## Journal de vérification *(2026-06-10 — ce qui a changé vs le premier jet)*

**Vérifs locales (Julia 1.12.5, chargement réel du module)** :
| Claim | Verdict | Preuve |
|---|---|---|
| `vcov_cluster`/`ClusterResult` exportés mais non définis | **CONFIRMÉ (bug)** | `isdefined(...)==false` → §0.1 |
| `m.sqrtwts == √w` (pondéré) / vide (non pondéré) | **CONFIRMÉ** | scaling ch.1 valide |
| `m.optsum.returnvalue` accessible (statut convergence) | **CONFIRMÉ** | `XTOL_REACHED` |
| `response(m)`=y, `scaled_re_matrix` renvoie ZΛ brut | **CONFIRMÉ** | ch.1 math |
| Manifest périmé après ajout `Random` | **CONFIRMÉ** | warning Pkg → §0.3 |

**Vérifs externes (sources citées)** :
| Claim (premier jet) | Verdict | Correction |
|---|---|---|
| clubSandwich types CR0/CR1/CR2 | **PARTIEL** | types réels : CR0, CR1, CR1p, **CR1S**, CR2, CR3 ; clubSandwich gère `lmer` ; **cluster-only, non spatial** |
| « CR1 = (G/(G−1))·((N−1)/(N−K)) » | **FAUX** | c'est **CR1S** (défaut Stata) ; CR1 = `m/(m−1)` seul → ch.2 corrigé |
| « Bartlett PSD-admissible (KP2007) » (package livré) | **SURÉVALUÉ** | `K₁∈P₁` (PSD 1-D only) ; 2-D exige `K₂∈P₂` → §0.2 + ch.3 |
| « uniforme non-PSD, *unlike Bartlett* » | **TROMPEUR** | en 2-D ni l'un ni l'autre garanti ; ancrer sur Schoenberg/Bochner |
| Référence « Huang 2022 » (CR0 mixte) | **MAL ATTRIBUÉ** | remplacer par **Huang, Wiedermann & Zhang 2022, MBR** et/ou Liang-Zeger 1986 |
| `plm::vcovSCC` = Driscoll-Kraay 1998 | **CONFIRMÉ** | ancre valide ch.6 ; DK1998 REStat 80(4):549-560 ✓ |
| CGM multiway `V=V₁+V₂−V₁₂` | **CONFIRMÉ** | ch.2 multiway OK |
| γ = semivariogramme = ½E[Δ²] ; variogramme strict = 2γ | **CONFIRMÉ** | ch.3.5 : renvoyer γ, pas 2γ |
| `gstat`/`GeoStats.jl` renvoient γ (semivariance), pas 2γ | **CONFIRMÉ** | pas de facteur 2 à la comparaison |
| `γ(h)=Ĉ(0)−Ĉ(h)` (sous sill fini) | **CONFIRMÉ** | property-test ch.3.5 + garde-plateau |
| « covariogramme » = C(h), usage Lehner cohérent | **PARTIEL** | terme surchargé (sens géométrie existe) ; OK en contexte |
| `gstat::variogram(z~1,d)` → np/dist/gamma ; `covariogram=TRUE`→C(h) | **CONFIRMÉ** | ⚠️ colonne reste `gamma` en mode covariogramme |
| `GeoStats.jl EmpiricalVariogram(...; nlags, maxlag, estimator)`, `values(γ)` | **CONFIRMÉ** | ⚠️ `nbins`→`nlags` (v0.7.0) ; **pas** de covariance empirique → Ĉ(h) via gstat seulement |
| JOSS paper « 250-1000 mots » | **FAUX** | **750-1750 mots** ; sections étendues + divulgation IA |
| JOSS seuil dur « 1000 lignes » | **OBSOLÈTE** | plus de gate LOC ; critères qualitatifs (effort/feature-complete/≥6 mois/citable) |
| `CoefTable` dans StatsAPI, champs nommés | **FAUX** | dans **StatsBase**, constructeur **positionnel** → ch.5 corrigé |
| Enregistrement General Registry = gate JOSS | **PARTIEL** | norme attendue, mais gate dur = release taggée + DOI Zenodo |

**Sources** : CRAN clubSandwich v0.7.0 & plm v2.6-7 ; Kelejian-Prucha 2007 (J.Econometrics 140, Remark 1) ; Cameron-Miller 2015 (JHR) ; Huang-Wiedermann-Zhang 2022 (MBR, doi:10.1080/00273171.2022.2077290) ; Driscoll-Kraay 1998 (REStat) ; joss.readthedocs.io (paper.md, review_criteria, rév. 2026) ; StatsAPI.jl & StatsBase.jl source ; gstat (r-spatial.github.io/gstat) ; GeoStatsFunctions.jl source ; Cressie 1993 (convention γ/2γ) ; Wikipedia Variogram (nugget/sill/range, γ=C(0)−C(h)).
