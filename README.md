# global_plant_ISAR

Code and data for the paper: The global island species–area relationship for plants
Authors: Thomas J. Matthews, Julian Schrader, François Rigal, Kostas A. Triantis, Holger Kreft, Patrick Weigelt, Robert J. Whittaker


## Datasets meta-data

Data were derived from Schrader et al. (2024), which is based on the Global Inventory of Flora and Traits (GIFT) database. Native species richness reflects the total number of species compiled from all available sources for each island.	

### data_global_SAR
island_name	Name of islands as given in GIFT
entity_class	binary indicating whether geographic entity is an island or continent
area	islan darea in km2
native_count	 richness of species native to the archipelago
endemic_count	number of species endemic to the island (Single Island Endemics, SIE)
category	oceanic (i.e. volcanic oceanic islands and atolls), continental (continental shelf), continental fragment
dist	distance of islabd from nearest mainland in km
longitude	centroid longitude of entity
latitude	centroid latitude of entity

### data_archipelago
archipelago_name	name of archipelgo as supplied in Schrader et al 2024
area	area of archipelago in km2
dist	distance of islands closest to the nearest mainland (km) in each archipelago
archipelago_cat	Predominant island categories within each archipelago: oceanic (i.e., volcanic oceanic archipelagos and atolls), continental (continental shelf islands), and continental fragments.
latitude	centroid longitude of archipelago
longitude	centroid latitude of archipelago
native_count	 species richness of species native to the archipelago
endemic_count	species richness of species endemic to the archipelago
![image](https://github.com/user-attachments/assets/d383c1de-d615-4924-b222-c540f3ea8c40)





