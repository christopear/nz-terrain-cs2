# Takaka river route snapshot

Retrieved 26 October 2026. Coordinates are longitude, latitude (WGS84).

© OpenStreetMap contributors, via WaterwayMap.org. The route database and
adaptations are available under the Open Database Licence (ODbL) 1.0:
https://www.openstreetmap.org/copyright
https://opendatacommons.org/licenses/odbl/1-0/

Download sources:

- Takaka: https://waterwaymap.org/river/T%C4%81kaka%20River%20000060780041/geometry.geojson
- Waingaro: https://waterwaymap.org/river/Waingaro%20River%20000060780046/geometry.geojson
- Anatoki: https://waterwaymap.org/river/Anatoki%20River%20000060780041/geometry.geojson
- Waikoropupu: https://waterwaymap.org/river/Te%20Waikoropup%C5%AB%20River%20001156528462/geometry.geojson

`takaka-routes.geojson` retains the longest connected LineString for each named
river. A short disconnected upper Anatoki fragment was omitted. `rivers.py`
clips routes to the playable footprint, snaps tributary outlets to the main
river only when within 5 m, and applies the explicit Waikoropupu upstream trim.
Route direction is upstream to downstream for this snapshot, checked against
terrain profiles and the connected river mouths. These are mapped centrelines,
not surveyed channel boundaries. Subsequent route edits require visual review.

Useful primary context:

- Tasman District Council, Takaka freshwater management unit:
  https://www.tasman.govt.nz/my-region/environment/environmental-management/water/water-management/freshwater-management-units-2/takaka-fmu
  Confirms the tributary network and Springs River's contribution to Waikoropupu.
- Council Takaka River Modelling Report, April 2024:
  https://www.tasman.govt.nz/document/serve/Takaka%20River%20Modelling%20Report%20April%202024.pdf?DocID=36203
  Describes supplementing LiDAR with 2022 Takaka/Waingaro bathymetry and 2019
  Anatoki cross-sections. Those survey datasets are NOT included in this project.

Synthetic channel widths and depths are experiment settings. They were not
inferred as measured hydraulic dimensions from those reports.

## Jamestown network

`jamestown-network.geojson` combines public LINZ river name line features with
complete main-river OSM centrelines. Snapshot downloaded 26 October 2026.

- LINZ River Name Lines Pilot, CC BY 4.0:
  https://services.arcgis.com/xdsHIIxuCWByZiCB/ArcGIS/rest/services/LINZ_NZ_River_Name_Lines_Pilot/FeatureServer/0
- Hollyford, OSM/WaterwayMap, ODbL 1.0:
  https://waterwaymap.org/river/Hollyford%20River%20/%20Whakatipu%20K%C4%81%20Tuka%20000060773327/geometry.geojson
- Awarua, OSM/WaterwayMap, ODbL 1.0:
  https://waterwaymap.org/river/Awarua%20River%20001050002387/geometry.geojson
- Dry Awarua, OSM/WaterwayMap, ODbL 1.0:
  https://waterwaymap.org/river/Dry%20Awarua%20River%20011933740643/geometry.geojson

Named LINZ Awarua/Dry Awarua fragments are replaced by complete OSM routes.
The combined route database and adaptations are available under ODbL 1.0;
retain both LINZ and OpenStreetMap attribution. Trust OSM downstream coordinate
order for these three added main routes. Other flow directions are inferred.
The source does not provide measured river widths: selection widths are an
explicit upstream-channel-length heuristic, with documented lake/border and
named assumptions. A five-metre threshold is not a field measurement.

Lake boundaries: LINZ Topo50 Lake Polygons, CC BY 4.0:
https://services.arcgis.com/xdsHIIxuCWByZiCB/ArcGIS/rest/services/NZ_Lake_Polygons_Topo_1_50k/FeatureServer/0
Lagoon boundaries: publicly hosted LINZ-derived Topo50 lagoon polygons (older
third-party distribution; see per-build source snapshot):
https://services6.arcgis.com/ZVM1rEuVZjtC1Wwk/ArcGIS/rest/services/Hydrography_Data/FeatureServer/18

All sea, lake, lagoon and river beds are synthetic. The McKerrow depth parameter
comes from Wellman and Wilson (1965), "Search for isostatically trapped sea water
in New Zealand": https://www.tandfonline.com/doi/pdf/10.1080/00288306.1965.10426432
This is a depth constraint only, not a source of surveyed basin geometry.


## Additional map networks (26 October 2026)

Port Underwood uses the public LINZ River Name Lines Pilot.
Milford Sound, Jackson Bay, Okarito and Bay of Islands combine that layer with
OSM/WaterwayMap main channels. Each added geometry retains its source URL.
Tūtoko/Tutoko aliases are reconciled to avoid duplicate centreline carving.
Main channels are synthetic narrow playable channels, not the full width of
braided river floodplains. No surveyed widths or discharge are claimed.

The centreline snapshots are inputs, not a guarantee of a complete or current
network. The exported selection report records trimming and exclusions.
