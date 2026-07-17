# Robotica Labo — studentenrepository

MATLAB-project voor de labozittingen robotica. Open het project via
`robotica_labo_student.prj`; dat stelt automatisch de juiste paden in.

## Aan de slag

1. Kloon de repository met **GitHub Desktop**:
   - Installeer [GitHub Desktop](https://desktop.github.com/).
   - Klik op deze GitHub-pagina op de groene knop **Code** en kies
     **Open with GitHub Desktop**, en bevestig de kloonlocatie.
   - GitHub Desktop haalt de submodules automatisch mee binnen.

   *Alternatief via de command line* (inclusief submodules):

   ```
   git clone --recurse-submodules https://github.com/sadegroo/robotica_labo_student.git
   ```

   (Al gekloond zonder submodules? Voer dan `git submodule update --init` uit.)

2. Begin met **`Voorbeelden/Getting_started.mlx`** — dit live script loodst je
   door de basis.

3. In **`Voorbeelden/`** vind je uitgewerkte voorbeelden (o.a.
   `MRvoorbeeld_4_1.m`) die je als vertrekpunt voor de oefeningen kunt
   gebruiken. De oefeningen zelf staan in `Oefeningen/`.

## Inbegrepen bibliotheken (bronnen)

| Map | Bibliotheek | Bron |
|---|---|---|
| `petercorke/rtb` | Robotics Toolbox for MATLAB (Peter Corke) | https://github.com/petercorke/robotics-toolbox-matlab |
| `petercorke/smtb` | Spatial Math Toolbox (Peter Corke) | https://github.com/petercorke/spatial-math |
| `petercorke/common` | Gemeenschappelijke hulpfuncties (Peter Corke) | https://github.com/petercorke/toolbox-common-matlab |
| `modernrobotics` | Modern Robotics-bibliotheek (Lynch & Park) | https://github.com/sadegroo/ModernRobotics (fork van https://github.com/NxRLab/ModernRobotics) |
| `MATLAB_Functions/arrow3D_pub` | arrow3D (Shawn Arseneau) | MATLAB File Exchange |
| `MATLAB_Functions/ST_pub` | Structure-tensor-demo | MATLAB File Exchange |
| `MATLAB_Functions/hypersphere` | hypersphere (Michael Völker) | MATLAB File Exchange (zie `license.txt`) |

### Documentatie Modern Robotics-functies

Een beschrijving van alle functies van de Modern Robotics-bibliotheek
(installatie, gebruik en de functielijst per hoofdstuk van het boek) vind je in
**`modernrobotics/packages/MATLAB/README.md`**. Na het openen van het project
kun je ook `help <functienaam>` gebruiken (bv. `help RotInv`).

## Eigen hulpfuncties

De overige functies in `MATLAB_Functions/` en `Extra_functions/` zijn eigen
hulpfuncties voor dit labo (o.a. `mat2latex`, `plotFrame`). Overzicht van
`Extra_functions/` (details: `help <functienaam>`):

**Algemeen**

| Interface | Doel |
|---|---|
| `[R, d] = nearestSO3(M)` | dichtstbijzijnde rotatiematrix + afstand |
| `[T, d] = nearestSE3(Ttilde)` | dichtstbijzijnde homogene transformatie + afstand |
| `e = manipulabilityEllipsoid(J, part)` | ellipsoïdedata (assen, stralen, w, conditie) uit een 6xn (`part` = `'linear'`/`'angular'`), 3xn of 2xn Jacobiaan |
| `h = plotEllipsoid(e, center, ...)` | tekent die ellipsoïde (3D-oppervlak of 2D-ellips) op positie `center` |
| `[e, fig] = ellips2R(theta, L)` | plot een planaire 2R arm (2 lijnsegmenten) met de manipuleerbaarheidsellips op de eindeffector en geeft de manipuleerbaarheidsmaten terug |

**`fwdkin/` — voorwaartse kinematica**

| Interface | Doel |
|---|---|
| `Rx = RotMatXaxis(angle, outputsize)` | elementaire rotatie om x (idem `RotMatYaxis`, `RotMatZaxis`) |
| `R = RotMatKaxis(k, theta)` | rotatie om willekeurige as k |
| `R = Euler2RotMat(angles, convention)` | rotatiematrix uit Euler-/fixed-angles (24 conventies, ook symbolisch) |
| `T = DH(alpha, a, d, theta)` | één klassieke DH-transformatie |
| `[Tfull, Tparts, Tcumul] = DH_full(DH_table)` | volledige keten, klassieke DH (idem `DH_modified`) |
| `T = fk_classicalDH(DH_table, theta)` | voorwaartse kinematica uit DH-tabel |
| `[Jv, Jw] = geometric_jacobian(DH_table, joint_types, frame_idx)` | geometrische Jacobiaan |
| `[linvel, angvel, Jv, Jw, symbols] = prop_vel_jac(DH_table, joint_types, ...)` | snelheidspropagatie |
| `D = translate(v)`, `Tinv = inv4(T)` | translatiematrix; snelle inverse van T |

**`invkin/` — inverse kinematica**

| Interface | Doel |
|---|---|
| `[sol1, sol2] = RotMat2Euler(R, convention)` | Euler-hoeken uit R (inverse van `Euler2RotMat`) |
| `[K1, th1, K2, th2] = RotMat2AxisAngle(R)` | as-hoekvoorstelling uit R |
| `solutions = pieperIK_classicalDH(DH_table, T_desired)` | Pieper-IK, klassieke DH (idem `pieperIK_modifiedDH`) |

**`tests/`** — unittests; uitvoeren met `runtests('Extra_functions/tests')`.
