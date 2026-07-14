# Robotica Labo — studentenrepository

MATLAB-project voor de labozittingen robotica. Open het project via
`robotica_labo_student.prj`; dat stelt automatisch de juiste paden in.

## Aan de slag

1. Kloon de repository **inclusief submodules**:

   ```
   git clone --recurse-submodules https://github.com/sadegroo/robotica_labo_student.git
   ```

   (Al gekloond zonder submodules? Voer dan `git submodule update --init` uit.)

2. Begin met **`Oefeningen/Getting_started.mlx`** — dit live script loodst je
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

De overige functies in `MATLAB_Functions/` en `Extra_functions/` zijn eigen
hulpfuncties voor dit labo (o.a. `mat2latex`, `plotFrame`, `nearestSO3`,
`nearestSE3`).
