# System Template

This is a template directory for an example system in the BEST-CSP project.
There should be a csv file for each polymorph, and each physical property for that polymorph, containing the experimental results,
alongside appropriate metadata for how the data was generated.
There is absolutely no expectation that every property will be filled in for each system,
but please try to include as much detail as possible for the data that is recorded.
The first line in the csv will be the REFCODE of the best representative structure for that system. The
second line will be the column titles, typically taking the format:
Identifier, Property, Value, Std, N, Name, Comment, Reference

- Identifer: A persistent identifier so that users can discuss individual results
- Property: The property, e.g. Melting point
- Value (Units): The measurand recorded in the experiment. Units will typically be K, kJ/mol
- Std: The Standard deviation of the value, if the contribution is a full dataset
- N: The number of experiments performed if the contribution is a full dataset
- Name: Contributor's name
- Reference/Author: If the value is from the literature, then leave a DOI, otherwise, an identifier for the set of unpublished data.
For example, I have used M_Dudek_2024.1 for the data contributed by Marta Dudek at the Zagreb meeting in March 2024
There should be a README.md in each System folder, with high level details about the system
- Comments: Anything the contributor considers pertinent to the reported value (pressure, ambient conditions etc.); please avoid commas and quote marks!

## Submitting data

In an ideal world, data would be contributed by trackers, branching from this github repository,
adding the data to the relevant place, and then submitting a Pull review, setting Isaac Sugden, Jonas Nyman,
Doris Braun, Ivo Rietveld, and the working group lead as reviewers. However, it is fine to submit data via email
and we can update the repo, but please try to include metadata either way

## Naming Conventions

csv files should be named in a way that unequivocably identifies the polymorph, e.g., Monoclinic_beta_II_melting_point.csv
if there is more than one way of referring to the polymorph

## README.md contents

within each README.md, name of compound, SMILES, pictures etc. could be recorded
