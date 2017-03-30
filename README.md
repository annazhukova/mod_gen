# mod_gen: Knowledge-based generalization for metabolic models

**Model Generalization** is a Python library that compresses a metabolic network model
using the knowledge-based model generalization method.

**Model generalization** takes a model in [SBML format](http://sbml.org/) as input, and produces 2 SBML files as an output:
* SBML containing the generalized model
* SBML file with [groups extension](http://sbml.org/Documents/Specifications/SBML_Level_3/Packages/groups)
containing the initial model plus the groups representing similar metabolites and similar reactions.

## Article

Zhukova A, Sherman DJ. **Knowledge-based generalization of metabolic models.**
*J Comput Biol.* 2014 Jul; **21**(7):534-47 [doi:10.1089/cmb.2013.0143](http://identifiers.org/doi/10.1089/cmb.2013.0143)


## Model Generalization Method

The **model generalization method** groups similar metabolites and reactions in the network
based on its structure and the knowledge, extracted from [the ChEBI ontology](http://www.ebi.ac.uk/chebi/).
The reactions between the same generalized species are factored together into generalized reactions.
A generalization is made specifically for a given model and is maximal with respect to the relations in the model;
it respects semantic constraints such as reaction stoichiometry, connectivity, and transport between compartments;
and it is performed through a heuristic method that is efficient in practice for genome-scale models.

Each metabolite can be generalized up to one of its ancestors in ChEBI. If a ChEBI annotation for a metabolite
is not present in the model, the method attempts to automatically deduce it by comparing metabolite’s name
to ChEBI terms’ names and synonyms. Reactions that share the same generalized reactants and
the same generalized products, are considered equivalent and are factored together into a generalized reaction.

The appropriate level of abstraction for metabolites and reactions is defined by the network itself as
the most general one that satisfies two restrictions:

1. *Stoichiometry preserving restriction:* metabolites that participate in the same reaction cannot be grouped together;

2. *Metabolite diversity restriction:* metabolites that do not participate in any pair of similar reactions are not
  grouped together (as there is no evidence of their similarity in the network).

Overall, the generalization method is composed of three modules:

1. Aggressive reaction grouping based on the most general metabolite grouping (defined by ChEBI),
in order to generate reaction grouping candidates;

2. Ungrouping of some metabolites and reactions to correct for violation of the stoichiometry preserving restriction;

3. Ungrouping of some metabolites (while keeping the reaction grouping intact) to correct for violation of
the metabolite diversity restriction.

For instance, *(S)-3-hydroxydecanoyl-CoA*, *(S)-3-hydroxylauroyl-CoA* and *(S)-3-hydroxytetradecanoyl-CoA*
have a common ancestor *hydroxy fatty acyl-CoA* in ChEBI. They can be grouped and generalized into *hydroxy fatty acyl-CoA*,
if in the network there is no reaction whose stoichiometry would be changed by such a generalization
(stoichiometry preserving restriction), and exist similar reactions that consume or produce them
(metabolite diversity restriction).


## Installation

From the directory where you have extracted this archive, execute:
```bash
python setup.py
```

## Running Model Generalization

Execute:

```bash
python3 ./sbml_generalization/runner/main.py --model path_to_your_model.xml --verbose
```

The script will produce two SBML files, containing the generalized model:

* path_to_your_model_generalized.xml -- SBML containing the generalized model
* path_to_your_model_with_groups.xml -- SBML file with groups extension containing the initial model
  plus the groups representing similar metabolites and similar reactions.
