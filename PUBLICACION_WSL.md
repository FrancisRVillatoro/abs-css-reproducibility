# Publicación desde WSL — GitHub y Zenodo

Esta guía presupone que el ZIP **GitHub** ha sido extraído en
`/mnt/c/Users/FRVillatoro/Downloads/abs-css-reproducibility`.
Los comandos NO publican nada automáticamente hasta ejecutar `git push`.

## 1. Verificación local antes de publicar

```bash
cd /mnt/c/Users/FRVillatoro/Downloads/abs-css-reproducibility
python3 scripts/audit_package.py
python3 repro/code/check_manifest.py
```

Si la ruta real de Downloads es distinta, corrige únicamente el primer `cd`.
No hace falta activar el entorno científico para estas comprobaciones.

## 2. Crea en GitHub un repositorio vacío

Nombre sugerido: `abs-css-reproducibility`. Propietario previsto:
`FrancisRVillatoro`. Al crear el repositorio, **NO** marques las casillas de
README, licencia, ni `.gitignore`: todos esos archivos ya existen aquí.
Elige visibilidad pública solo cuando hayas revisado el contenido y la licencia.

## 3. Git desde WSL

```bash
git init -b main
git add .
git status --short
git commit -m "Release reproducibility code and data v1.0.0"
git remote add origin https://github.com/FrancisRVillatoro/abs-css-reproducibility.git
git push -u origin main
```

Si `git remote add origin` indica que el remoto ya existe, **no repitas el
comando a ciegas**. Comprueba `git remote -v` y actualiza solo si corresponde.

## 4. GitHub → Zenodo (método recomendado)

1. Vincula tu cuenta de GitHub en Zenodo (perfil → Linked accounts).
2. En Zenodo, perfil → GitHub → **Sync now** y activa este repositorio.
3. **Después de activarlo**, crea el tag y súbelo desde WSL:

```bash
git tag -a v1.0.0 -m "v1.0.0 — audited code and data"
git push origin v1.0.0
```

4. En GitHub, abre **Releases → Draft a new release**, selecciona `v1.0.0`
   y publica la release.
5. Espera a que Zenodo procese el release y comprueba en Zenodo el título,
   autor, versión, licencia, archivos y DOI definitivo.

**Importante:** un `git push` del tag no equivale a publicar una GitHub Release.
Zenodo suele archivar la release, no un mero commit/tag.

## 5. Alternativa: depósito manual Zenodo

Si NO usarás la integración GitHub→Zenodo, puedes cargar directamente
`abs-css-reproducibility-v1.0.0-zenodo.zip` en **New upload**.
Recurso: **Software**; autor: Villatoro, Francisco R.; versión: 1.0.0;
acceso Open; licencia MIT (revisar derechos de distribución);
fecha 2026-10-10. Publica solo cuando los metadatos estén verificados.

**No hagas integración automática y depósito manual del mismo release** salvo
que quieras generar deliberadamente registros/DOI duplicados.

## 6. Verificación final

El repositorio GitHub y la versión Zenodo deben incluir código, CSV/JSON/NPZ,
certificados y metadatos, **no manuscrito, paper, LaTeX, PDF, figuras PDF/PNG
ni carta editorial**. El script `scripts/audit_package.py` verifica esta regla.

Cuando Zenodo asigne el DOI, consérvalo para la Data Availability Statement del
paper. El paquete no contiene ese manuscrito y no intenta editarlo.
