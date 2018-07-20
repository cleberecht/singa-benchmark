package bio.singa.structure;

import bio.singa.core.utility.Resources;
import bio.singa.structure.parser.pdb.structures.StructureParserOptions;

import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * @author fk
 */
public interface BenchmarkConstants {
    String LOCAL_PDB_LOCATION = "/srv/pdb";
    Path CHAIN_LIST_PATH_500 = Paths.get(Resources.getResourceAsFileLocation("nrpdb_041416_BLAST_e-7_pdb-mmtf_mutual_subset_500.txt"));
    Path CHAIN_LIST_PATH_100 = Paths.get(Resources.getResourceAsFileLocation("nrpdb_041416_BLAST_e-7_pdb-mmtf_mutual_subset_100.txt"));
    Path CHAIN_LIST_PATH_10 = Paths.get(Resources.getResourceAsFileLocation("nrpdb_041416_BLAST_e-7_pdb-mmtf_mutual_subset_10.txt"));
    StructureParserOptions STRUCTURE_PARSER_OPTIONS = StructureParserOptions.withSettings(StructureParserOptions.Setting.OMIT_EDGES,
                                                                                          StructureParserOptions.Setting.OMIT_HYDROGENS,
                                                                                          StructureParserOptions.Setting.OMIT_LIGAND_INFORMATION,
                                                                                          StructureParserOptions.Setting.GET_IDENTIFIER_FROM_FILENAME);
}
