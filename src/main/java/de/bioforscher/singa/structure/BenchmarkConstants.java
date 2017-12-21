package de.bioforscher.singa.structure;

import de.bioforscher.singa.core.utility.Resources;

import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * @author fk
 */
public interface BenchmarkConstants {
    String LOCAL_PDB_LOCATION = "/srv/pdb";
    Path CHAIN_LIST_PATH = Paths.get(Resources.getResourceAsFileLocation("nrpdb_041416_BLAST_e-7_pdb-mmtf_mutual_subset_500.txt"));
}
