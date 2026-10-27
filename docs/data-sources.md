# Data sources and acquisition

All source downloads are public and require no paid bathymetry or account keys.
Availability and mapping completeness vary by location. The supplied snapshots
were assembled for the configured presets; check each package's manifest and
water/river reports for the sources actually used.

## Elevations

| Source | Role | Attribution / reference |
|---|---|---|
| LINZ public LiDAR DEM collections | Preferred bare-earth elevations; selected explicitly per preset, later collections take precedence | [STAC catalogue](https://nz-elevation.s3.ap-southeast-2.amazonaws.com/catalog.json), [archive documentation](https://github.com/linz/elevation); LINZ/provider attribution and CC BY 4.0 as documented by each collection |
| LINZ NZ 8 m DEM | Fallback terrain interpolated from Topo50 contours; 8 m cells do not imply 8 m surveyed detail | [Dataset metadata](https://data.linz.govt.nz/layer/51768-nz-8m-digital-elevation-model-2012/); Geographx / LINZ, CC BY 4.0 |
| Copernicus GLO-30 Public | Surface-elevation fallback; contains vegetation/buildings | [Dataset and licence links](https://registry.opendata.aws/copernicus-dem/), [format documentation](https://copernicus-dem-30m.s3.amazonaws.com/readme.html); retain the Copernicus/DLR/Airbus acknowledgement in generated packages |

The builder queries intersecting source tiles and reads raster windows. Source
URLs, collection information and output checksums are retained in manifests.
Copernicus ocean tiles can be absent: a 404 is filled with zero only after the
tile is absent from both GLO-30 and worldwide GLO-90 catalogues. Other failures
remain errors. No paid bathymetry is downloaded.

## Water boundaries and river routes

- [LINZ Topo50 lake polygons, public ArcGIS distribution](https://services.arcgis.com/xdsHIIxuCWByZiCB/ArcGIS/rest/services/NZ_Lake_Polygons_Topo_1_50k/FeatureServer/0):
  lake boundaries and available surface levels, CC BY 4.0.
- [LINZ-derived lagoon polygons](https://services6.arcgis.com/ZVM1rEuVZjtC1Wwk/ArcGIS/rest/services/Hydrography_Data/FeatureServer/18):
  an older third-party distribution; retained source attribution and snapshots
  make this distinction explicit.
- [LINZ River Name Lines Pilot](https://services.arcgis.com/xdsHIIxuCWByZiCB/ArcGIS/rest/services/LINZ_NZ_River_Name_Lines_Pilot/FeatureServer/0):
  mapped river/stream centrelines, CC BY 4.0. Wide river channels can be missing
  from this line layer, and offshore-island coverage is incomplete.
- [OpenStreetMap](https://www.openstreetmap.org/copyright), via
  [WaterwayMap](https://waterwaymap.org/): supplementary main
  river channels. © OpenStreetMap contributors,
  ODbL 1.0; retain attribution and applicable database share-alike obligations.
  Adapted route/polygon geometry is supplied with packages.

See [the curated route source inventory](../river-configs/SOURCES.md) for
individual main-river downloads and references. Feature-level source URLs are
also stored in supplemental GeoJSON files. Do not infer measured widths, depths
or discharges from these centreline datasets.

Lake McKerrow's approximate 109 m depth parameter is sourced from
[Wellman and Wilson (1965)](https://www.tandfonline.com/doi/pdf/10.1080/00288306.1965.10426432).
Its basin shape is still synthetic. Defaults are 30 m for lakes and 10 m for
lagoons; narrow basins may never reach those ceilings.

## Preparing and rebuilding inputs

Existing curated configurations and route snapshots are already in the checkout.
Ordinary builds download missing elevations and water polygons automatically.
For deliberate input preparation:

```sh
uv run nz-terrain-discover-rivers --location jamestown
uv run nz-terrain --location jamestown --with-rivers
```

The LINZ discovery command queries the playable footprint and preserves existing
OSM supplements. It does not discover all missing wide channels automatically.
Takaka uses its existing four tuned
routes rather than the generic discovered network.

Acquisition commands reuse cached responses. To obtain fresh upstream data,
remove the relevant metadata/download cache deliberately, then prepare the inputs
again. Location changes generate different grid-cache keys. Rebuild both baseline
and river candidate after changing inputs; the original uncarved surface is
required so river depths are not stacked onto already carved lake beds.

`data/` and `output/` are excluded from Git. Small curated GeoJSON files in
`river-configs/` and `water-configs/`, plus test fixtures, remain tracked. Keep
their attribution when sharing or replacing them.
