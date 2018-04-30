package de.bioforscher.singa.structure.algorithms.superimposition.fit3d;

import de.bioforscher.singa.core.utility.Resources;
import de.bioforscher.singa.structure.BenchmarkConstants;
import de.bioforscher.singa.structure.model.oak.StructuralEntityFilter;
import de.bioforscher.singa.structure.model.oak.StructuralMotif;
import de.bioforscher.singa.structure.model.oak.Structures;
import de.bioforscher.singa.structure.parser.pdb.structures.SourceLocation;
import de.bioforscher.singa.structure.parser.pdb.structures.StructureParser;
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
import java.nio.file.Paths;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

/**
 * A benchmark case to test the runtime of the Fit3D algorithm for query motifs of extents. Thee extents are binned according to Math.floor(). The query motifs were selected from a non-redundant
 * snapshot (BLAST_e-80) of the Catalytic Site Atlas according to nrpdb.041416 VAST and have to meet the following criteria:
 * <pre>
 *     #1. motif is not allowed to contain repetitive amino acids
 *     2. maximal motif extent has to be <=20 Angstroem
 *     3. sorted by ascending extent only the first EXTENT_RANK_CUTOFF motifs are considered for the benchmark
 * </pre>
 *
 * @author fk
 */
@State(Scope.Benchmark)
public class Fit3DBenchmarkMotifByExtent {

    public static final int EXTENT_RANK_CUTOFF = 3;
    @Param({"5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"})
    private int motifExtent;
    private List<StructuralMotif> queryMotifs;

    public static void main(String[] args) throws RunnerException, IOException {
        Options opt = new OptionsBuilder()
                .include(Fit3DBenchmarkMotifByExtent.class.getSimpleName())
                .warmupIterations(5)
                .measurementIterations(5)
                .forks(1)
                .mode(Mode.AverageTime)
                .timeout(TimeValue.hours(24))
                .timeUnit(TimeUnit.MILLISECONDS)
                .build();
        Collection<RunResult> results = new Runner(opt).run();
        StringJoiner stringJoiner = new StringJoiner("\n", "extent,min,max,mean,stdev,ci95_min,ci95_max\n", "");
        for (RunResult result : results) {
            Statistics statistics = result.getPrimaryResult().getStatistics();
            double[] confidenceInterval = statistics.getConfidenceIntervalAt(0.95);
            String resultLine = result.getParams().getParam("motifExtent") +
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
        Files.write(Paths.get("results_extent.csv"), stringJoiner.toString().getBytes());
    }

    @Setup
    public void setUp() throws IOException {
        queryMotifs = Files.walk(Paths.get(Resources.getResourceAsFileLocation("structural_motifs/extent/extent_" + motifExtent)))
                           .filter(path -> path.toFile().isFile())
                           .map(path -> StructuralMotif.fromLeafSubstructures(StructureParser.local().path(path).parse().getAllLeafSubstructures()))
                           .sorted(Comparator.comparing(Structures::calculateExtent))
                           .limit(EXTENT_RANK_CUTOFF)
                           .collect(Collectors.toList());
    }

    @Benchmark
    @OperationsPerInvocation(EXTENT_RANK_CUTOFF)
    public void runWithLocalMmtf() {
        System.out.println("processing motif extent bin " + motifExtent);
        for (int i = 0; i < queryMotifs.size(); i++) {
            StructuralMotif queryMotif = queryMotifs.get(i);
            System.out.println("processing motif " + (i + 1) + "/" + queryMotifs.size());
            StructureParser.MultiParser multiParser = StructureParser.local()
                                                                     .localPDB(new StructureParser.LocalPDB(BenchmarkConstants.LOCAL_PDB_LOCATION, SourceLocation.OFFLINE_MMTF))
                                                                     .chainList(BenchmarkConstants.CHAIN_LIST_PATH_500, "\t")
                                                                     .setOptions(BenchmarkConstants.STRUCTURE_PARSER_OPTIONS);
            Fit3DBuilder.create()
                        .query(queryMotif)
                        .targets(multiParser)
                        .maximalParallelism()
                        .atomFilter(StructuralEntityFilter.AtomFilter.isArbitrary())
                        .run();
        }
    }
}
