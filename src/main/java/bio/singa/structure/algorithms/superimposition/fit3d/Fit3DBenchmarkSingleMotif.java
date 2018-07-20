package bio.singa.structure.algorithms.superimposition.fit3d;

import bio.singa.core.utility.Resources;
import bio.singa.structure.BenchmarkConstants;
import bio.singa.structure.model.families.AminoAcidFamily;
import bio.singa.structure.model.identifiers.LeafIdentifier;
import bio.singa.structure.model.oak.StructuralEntityFilter;
import bio.singa.structure.model.oak.StructuralMotif;
import bio.singa.structure.parser.pdb.structures.SourceLocation;
import bio.singa.structure.parser.pdb.structures.StructureParser;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.results.RunResult;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;
import org.openjdk.jmh.runner.options.TimeValue;
import org.openjdk.jmh.util.Statistics;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.StringJoiner;
import java.util.concurrent.TimeUnit;

/**
 * A benchmark case to test the runtime of the Fit3D algorithm for a single query motif, the catalytic triad of serine proteases. Benchmark is run for local PDB and MMTF installation, respectively.
 *
 * @author fk
 */
@State(Scope.Benchmark)
public class Fit3DBenchmarkSingleMotif {

    //    @Param({"PDB", "MMTF"})
    @Param({"MMTF"})
    private String parsing;

    @Param({"250", "500", "750", "1000", "1250", "1500", "1750", "2000"})
//    @Param({"250", "500"})
    private int datasetSize;
    private StructuralMotif queryMotif;

    public static void main(String[] args) throws RunnerException, IOException {
        Options opt = new OptionsBuilder()
                .include(Fit3DBenchmarkSingleMotif.class.getSimpleName())
                .warmupIterations(5)
                .measurementIterations(5)
                .forks(1)
                .mode(Mode.AverageTime)
                .timeout(TimeValue.hours(24))
                .timeUnit(TimeUnit.MILLISECONDS)
                .build();
        Collection<RunResult> results = new Runner(opt).run();
        StringJoiner stringJoiner = new StringJoiner("\n", "parsing,dataset_size,min,max,mean,stdev,ci95_min,ci95_max\n", "");
        for (RunResult result : results) {
            Statistics statistics = result.getPrimaryResult().getStatistics();
            double[] confidenceInterval = statistics.getConfidenceIntervalAt(0.95);
            String resultLine = result.getParams().getParam("parsing") +
                                "," +
                                result.getParams().getParam("datasetSize") +
                                "," +
                                statistics.getMin() +
                                "," +
                                statistics.getMax() +
                                "," +
                                statistics.getMean() +
                                "," +
                                statistics.getStandardDeviation() +
                                "," +
                                confidenceInterval[0] +
                                "," +
                                confidenceInterval[1];
            stringJoiner.add(resultLine);
        }
        Files.write(Paths.get("results_single_motif.csv"), stringJoiner.toString().getBytes());
    }

    @Setup
    public void setUp() {
//        queryMotif = StructuralMotif.fromLeafSubstructures(StructureParser.local()
//                                                                          .path(Paths.get(Resources.getResourceAsFileLocation("structural_motifs/1GL0_HDS_intra_E-H57_E-D102_E-S195.pdb")))
//                                                                          .parse()
//                                                                          .getAllLeafSubstructures());
        StructuralMotif structuralMotif = StructuralMotif.fromLeafSubstructures(StructureParser.local()
                                                                                               .inputStream(Resources.getResourceAsStream("structural_motifs/motif_KDEEH.pdb"))
                                                                                               .parse()
                                                                                               .getAllLeafSubstructures());
        structuralMotif.addExchangeableFamily(LeafIdentifier.fromSimpleString("A-164"), AminoAcidFamily.HISTIDINE);
        structuralMotif.addExchangeableFamily(LeafIdentifier.fromSimpleString("A-247"), AminoAcidFamily.ASPARTIC_ACID);
        structuralMotif.addExchangeableFamily(LeafIdentifier.fromSimpleString("A-247"), AminoAcidFamily.ASPARAGINE);
        structuralMotif.addExchangeableFamily(LeafIdentifier.fromSimpleString("A-297"), AminoAcidFamily.LYSINE);
        queryMotif = structuralMotif;
    }

    @Benchmark
    public void runForSingleMotif() {

        SourceLocation sourceLocation;
        if (parsing.equals("PDB")) {
            sourceLocation = SourceLocation.OFFLINE_PDB;
        } else {
            sourceLocation = SourceLocation.OFFLINE_MMTF;
        }

        String fileName = "nrpdb_041416_BLAST_e-7_pdb-mmtf_mutual_subset_" + datasetSize + ".txt";
        Path chainListPath = Paths.get(Resources.getResourceAsFileLocation(fileName));

        StructureParser.MultiParser multiParser = StructureParser.local()
                                                                 .localPDB(new StructureParser.LocalPDB(BenchmarkConstants.LOCAL_PDB_LOCATION, sourceLocation))
                                                                 .chainList(chainListPath, "\t")
                                                                 .setOptions(BenchmarkConstants.STRUCTURE_PARSER_OPTIONS);
        Fit3D run = Fit3DBuilder.create()
                                .query(queryMotif)
                                .targets(multiParser)
                                .maximalParallelism()
                                .atomFilter(StructuralEntityFilter.AtomFilter.isArbitrary())
                                .filterEnvironments(5)
                                .run();
        System.out.println(run.getMatches().size());
    }
}
