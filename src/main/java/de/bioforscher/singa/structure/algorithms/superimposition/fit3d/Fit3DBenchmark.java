package de.bioforscher.singa.structure.algorithms.superimposition.fit3d;

import de.bioforscher.singa.core.utility.Resources;
import de.bioforscher.singa.structure.BenchmarkConstants;
import de.bioforscher.singa.structure.model.identifiers.LeafIdentifiers;
import de.bioforscher.singa.structure.model.interfaces.Structure;
import de.bioforscher.singa.structure.model.oak.StructuralEntityFilter.AtomFilter;
import de.bioforscher.singa.structure.model.oak.StructuralMotif;
import de.bioforscher.singa.structure.parser.pdb.structures.SourceLocation;
import de.bioforscher.singa.structure.parser.pdb.structures.StructureParser;
import de.bioforscher.singa.structure.parser.pdb.structures.StructureParser.LocalPDB;
import de.bioforscher.singa.structure.parser.pdb.structures.StructureParser.MultiParser;
import de.bioforscher.singa.structure.parser.pdb.structures.StructureParserOptions;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;
import org.openjdk.jmh.runner.options.TimeValue;

import java.util.concurrent.TimeUnit;

/**
 * A benchmark case that is run to analyze the performance of the Fit3D algorithm.
 *
 * @author fk
 */
@State(Scope.Benchmark)
public class Fit3DBenchmark {

    private StructureParserOptions structureParserOptions;
    private StructuralMotif queryMotif;
    private Structure pdbTarget;
    private Structure mmtfTarget;

    public static void main(String[] args) throws RunnerException {
        Options opt = new OptionsBuilder()
                .include(Fit3DBenchmark.class.getSimpleName())
                .warmupIterations(5)
                .measurementIterations(5)
                .forks(1)
                .mode(Mode.AverageTime)
                .timeout(TimeValue.hours(1))
                .timeUnit(TimeUnit.MILLISECONDS)
                .build();
        new Runner(opt).run();
    }

    @Setup
    public void setUp() {
        Structure motifContainingStructure = StructureParser.local()
                                                            .fileLocation(Resources.getResourceAsFileLocation("1GL0_HDS_intra_E-H57_E-D102_E-S195.pdb"))
                                                            .parse();
        queryMotif = StructuralMotif.fromLeafIdentifiers(motifContainingStructure,
                                                         LeafIdentifiers.of("E-57", "E-102", "E-195"));

        pdbTarget = StructureParser.pdb()
                                   .pdbIdentifier("4cha")
                                   .parse();
        mmtfTarget = StructureParser.mmtf()
                                    .pdbIdentifier("4cha")
                                    .parse();

        structureParserOptions = StructureParserOptions.withSettings(StructureParserOptions.Setting.OMIT_EDGES,
                                                                     StructureParserOptions.Setting.OMIT_HYDROGENS,
                                                                     StructureParserOptions.Setting.OMIT_LIGAND_INFORMATION,
                                                                     StructureParserOptions.Setting.GET_IDENTIFIER_FROM_FILENAME);
    }

    @Benchmark
    public void runWithLocalPdb() {
        MultiParser multiParser = StructureParser.local()
                                                 .localPDB(new LocalPDB(BenchmarkConstants.LOCAL_PDB_LOCATION, SourceLocation.OFFLINE_PDB))
                                                 .chainList(BenchmarkConstants.CHAIN_LIST_PATH, "\t")
                                                 .setOptions(structureParserOptions);
        Fit3DBuilder.create()
                    .query(queryMotif)
                    .targets(multiParser)
                    .maximalParallelism()
                    .atomFilter(AtomFilter.isArbitrary())
                    .run();
    }

    @Benchmark
    public void runWithLocalMmtf() {
        MultiParser multiParser = StructureParser.local()
                                                 .localPDB(new LocalPDB(BenchmarkConstants.LOCAL_PDB_LOCATION, SourceLocation.OFFLINE_MMTF))
                                                 .chainList(BenchmarkConstants.CHAIN_LIST_PATH, "\t")
                                                 .setOptions(structureParserOptions);
        Fit3DBuilder.create()
                    .query(queryMotif)
                    .targets(multiParser)
                    .maximalParallelism()
                    .atomFilter(AtomFilter.isArbitrary())
                    .run();
    }

    @Benchmark
    public void runAgainstPdbStructure() {
        Fit3DBuilder.create()
                    .query(queryMotif)
                    .target(pdbTarget)
                    .atomFilter(AtomFilter.isArbitrary())
                    .run();
    }

    @Benchmark
    public void runAgainstMmtfStructure() {
        Fit3DBuilder.create()
                    .query(queryMotif)
                    .target(mmtfTarget)
                    .atomFilter(AtomFilter.isArbitrary())
                    .run();
    }
}