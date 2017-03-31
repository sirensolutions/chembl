pipeline changes
================
- switched from using posgres dumps + logstash to use sqlite dumps + python. The full process of importing a chembl relase is reproducible and takes now less an 1h.

data changes
============
- changed concat strings to be stored a lists. Many biological replicates has commas inside. this should improve the querying precision.
- added target type and uniprot accession to target index
- added assays type data
- fixed mappings for papers abstract and title, so they are mapped as english text and the word clouds show meaningful words.
- fixed mappings for papers year, not indexed as a date
- fixed mappings for assays description, so it is now mapped as english text and the word clouds show meaningful words.
- fixed mappings for molecule usan_year, not indexed as a date
- added compound records index
- added data to the molecule index (big extension)

kibi changes
============

- set replicas to 0 and shard to 1
- revised dashboards for:

* molecules
* targets
* papers
* assays
* activities


- Targets

  * improved table data
  * added external link to chembl
  * added external link to uniprot
  * added heatmap organism/target_type
  * added organism bubble chart

- Molecules

  * improved table data
  * added heatmap of molecule type to indication
  * added donut chart of phase to indication and source to indication


- Papers

  * improved table data
  * added cloud word for title
  * added cloud word for abstract

- Assays

  * improved table data
  * added percent barchart for assay score over assay type and tissue
  * added organism bubble diagram
  * added word cloud for assay description
  * added word cloud for assay type


- Activity

  * improved table data
  * added activity value box plot
