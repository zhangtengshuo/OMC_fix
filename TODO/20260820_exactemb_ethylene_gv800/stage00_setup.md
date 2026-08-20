# Stage 0: setup and provenance

Status: passed

Started: 2026-08-20T10:34:19+08:00

## Fixed experiment contract

- Run as `alpha` on host `GV800`.
- Write all new remote artifacts below `/scratch/alpha/exactemb-ethylene-20260820_103419`.
- Use `/scratch/alpha/OpenMolcasUP_ExactDenEmb-65a3eeb11` and preserve its dependent `/scratch/alpha/build-exactden-20260820_093910/wignernj-sys` installation.
- Use the fixed parallel ethylene dimer geometry with fragment A at z=0 Angstrom and fragment B at z=3.5 Angstrom.
- Use `cc-pVDZ`, C1, identical atom order, basis convention, and global orientation in fragment and dimer contexts.
- Use one fixed MPI process count for all RunFiles participating in an EXACTEMB construction.
- Keep the full-dimer Cholesky context separate from each target-fragment SEWARD/XField context.
- Treat partner nuclear attraction and partner electronic Coulomb embedding as separate records that must each appear exactly once.
- Never use an embedded Hamiltonian for final NOCI B3 integral generation.

## Initial observations

- `exactemb.exe` exists at `/scratch/alpha/OpenMolcasUP_ExactDenEmb-65a3eeb11/bin/exactemb.exe`.
- The installed `pymolcas` reports Python driver version `py2.32`.
- Slurm is available and the `molcas` partition has one 64-core node.
- The historical Dimerge ethylene results are principally below `/scratch/alpha/Dimerge/data/gv800_runs`.
- Historical results retain orbitals, VecDet files, HDF5 outputs, logs, and PySCF NPZ files, but no reusable RunFile, ONEINT, or Cholesky context was found.
- A cc-pVDZ ethylene monomer has 48 basis functions in the accepted Round 2 logs, so an A-then-B dimer is expected to have 96 AOs and candidate first-AO indices A=1 and B=49.
- The AO indices remain provisional until the regenerated dimer log and the EXACTEMB electron-count invariant confirm them.
- The historical `oneshot_R3p5_ccpvdz/oneshot.npz` contains isolated and embedded MO coefficients, embedded AO densities, and one-cycle history, but it does not contain the direct-J matrix.
- The older aug-cc-pVDZ ionic reference has a known unbound-anion artifact and will not be used as the first implementation gate.

## Stage 0 acceptance criteria

- Record exact host, user, OpenMolcas version/build, compiler/runtime libraries, scheduler, geometry, basis, process count, and remote directory layout.
- Establish a non-overwriting Slurm runner that keeps task-owned scratch below the experiment root and preserves requested RunFiles and integral contexts.
- Confirm the experiment begins with no running or stale job carrying the new experiment tag.

## Error and repair log

No experiment error has occurred.

## Provenance result

- Host and identity were confirmed as `GV800`, uid 1004 `alpha`.
- OpenMolcas reports `.molcasversion` value `o65a3eeb11dcbde ?`.
- `exactemb.exe` SHA-256 is `b4c8244f9d5b7d49588453c7af6644dfc642250cd4ad35ba0234b56e0b6747d3`.
- `pymolcas` SHA-256 is `eae30a6f23b7444fbac1b6ea2e1a44e986206a086837cb997ce7d10b1dd57712`.
- `ldd exactemb.exe` found `libmolcas`, HDF5, LibXC, MKL, and the offline wignernj libraries with no unresolved dependency.
- The compiler is GNU Fortran 13.3.0 and the MPI wrapper reports Open MPI 4.1.6.
- At 2026-08-20T10:36:10+08:00, no `alpha` job was running or queued and the `molcas` partition reported 64 idle CPUs.
- The remote directory layout and the local mirrored Markdown records were both created successfully.

Completed: 2026-08-20T10:36:10+08:00
