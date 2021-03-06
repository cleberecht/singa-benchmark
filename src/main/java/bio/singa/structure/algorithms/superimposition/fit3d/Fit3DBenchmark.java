package bio.singa.structure.algorithms.superimposition.fit3d;

import bio.singa.core.utility.Resources;
import bio.singa.structure.BenchmarkConstants;
import bio.singa.structure.model.identifiers.LeafIdentifiers;
import bio.singa.structure.model.interfaces.Structure;
import bio.singa.structure.model.oak.StructuralEntityFilter.AtomFilter;
import bio.singa.structure.model.oak.StructuralMotif;
import bio.singa.structure.parser.pdb.structures.SourceLocation;
import bio.singa.structure.parser.pdb.structures.StructureParser;
import bio.singa.structure.parser.pdb.structures.StructureParser.LocalPDB;
import bio.singa.structure.parser.pdb.structures.StructureParser.MultiParser;
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
                                                            .fileLocation(Resources.getResourceAsFileLocation("structural_motifs/1GL0_HDS_intra_E-H57_E-D102_E-S195.pdb"))
                                                            .parse();
        queryMotif = StructuralMotif.fromLeafIdentifiers(motifContainingStructure,
                                                         LeafIdentifiers.of("E-57", "E-102", "E-195"));

        pdbTarget = StructureParser.pdb()
                                   .pdbIdentifier("4cha")
                                   .parse();
        mmtfTarget = StructureParser.mmtf()
                                    .pdbIdentifier("4cha")
                                    .parse();
    }

    @Benchmark
    public void runWithLocalPdb() {
        MultiParser multiParser = StructureParser.local()
                                                 .localPDB(new LocalPDB(BenchmarkConstants.LOCAL_PDB_LOCATION, SourceLocation.OFFLINE_PDB))
                                                 .chainList(BenchmarkConstants.CHAIN_LIST_PATH_500, "\t")
                                                 .setOptions(BenchmarkConstants.STRUCTURE_PARSER_OPTIONS);
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
                                                 .chainList(BenchmarkConstants.CHAIN_LIST_PATH_500, "\t")
                                                 .setOptions(BenchmarkConstants.STRUCTURE_PARSER_OPTIONS);
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