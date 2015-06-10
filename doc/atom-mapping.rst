.. _atom-mapping:

Residue-residue to atom-atom mapping
====================================

The result of the DCA analysis is a list of nucleotide contacts. These need to
be mapped to atom-atom contacts to be usable as constraints in Rosetta. Two
options are implemented.

P-Only mapping: ``pOnly``
-------------------------

A simple variant is to take the P atom of each nucleotide backbone and add
these as atom-atom constraints.

Complex atom-atom mapping using a contact database: ``allAtomWesthof``
----------------------------------------------------------------------

This method relies on pre-built contact lists included in the prediction
package. For each type of nucleotide contact (A-A, A-C, A-G, A-U, ...) a list of
representatives is stored, containing the PDB code of an entry with known
native structure and the numbers of the residues. From these lists a database
is built containing all possible atom-to-atom combinations for the two residues,
as well as their distances extracted from the crystal structures. Next, to
reduce this database to the relevant atoms involved in forming a contact, a
second database is created, containing only those atom-to-atom contacts that
satisfy a mean and standard deviation condition. The defaults for those are set
to lower than 6Å and 3Å respectively. DCA contacts such as "A-U" are then
looked up in the mean database, and all atom-to-atom combinations found are
added as constraints.
