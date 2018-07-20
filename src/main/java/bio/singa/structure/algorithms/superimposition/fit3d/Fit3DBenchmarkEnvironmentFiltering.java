package bio.singa.structure.algorithms.superimposition.fit3d;

import bio.singa.core.utility.Resources;
import bio.singa.structure.BenchmarkConstants;
import bio.singa.structure.model.interfaces.LeafSubstructure;
import bio.singa.structure.model.oak.StructuralEntityFilter;
import bio.singa.structure.model.oak.StructuralMotif;
import bio.singa.structure.model.oak.Structures;
import bio.singa.structure.parser.pdb.structures.SourceLocation;
import bio.singa.structure.parser.pdb.structures.StructureParser;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.results.RunResult;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;
import org.openjdk.jmh.util.Statistics;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.StringJoiner;
import java.util.concurrent.TimeUnit;

/**
 * A benchmark case to test the runtime of the Fit3D algorithm for query motifs of different size. The query motifs were selected from a non-redundant snapshot (BLAST_e-80) of the
 * Catalytic Site Atlas according to nrpdb.041416 VAST and have to meet the following criteria:
 * <pre>
 *     #1. motif is not allowed to contain repetitive amino acids
 *     2. maximal motif extent has to be <=20 Angstroem
 *     3. motif is not allowed to contain more than 6 amino acids
 * </pre>
 *
 * @author fk
 */
@State(Scope.Benchmark)
public class Fit3DBenchmarkEnvironmentFiltering {

    @Param({
                   "286",
                   "285",
                   "284",
                   "283",
//            "282",
//            "281",
//            "280",
//            "279",
//            "278",
//            "277",
//            "276",
//            "275",
//            "274",
//            "273",
//            "272",
//            "271",
//            "270",
//            "269",
//            "268",
//            "267",
//            "266",
//            "265",
//            "264",
//            "263",
//            "262",
//            "261",
//            "260",
//            "259",
//            "258",
//            "257",
//            "256",
//            "255",
//            "254",
//            "253",
//            "252",
//            "251",
//            "250",
//            "249",
//            "248",
//            "247",
//            "246",
//            "245",
//            "244",
//            "243",
//            "242",
//            "241",
//            "240",
//            "239",
//            "238",
//            "237",
//            "236",
//            "235",
//            "234",
//            "233",
//            "232",
//            "231",
//            "230",
//            "229",
//            "228",
//            "227",
//            "226",
//            "225",
//            "224",
//            "223",
//            "222",
//            "221",
//            "220",
//            "219",
//            "218",
//            "217",
//            "216",
//            "215",
//            "214",
//            "213",
//            "212",
//            "211",
//            "210",
//            "209",
//            "208",
//            "207",
//            "206",
//            "205",
//            "204",
//            "203",
//            "202",
//            "201",
//            "200",
//            "199",
//            "198",
//            "197",
//            "196",
//            "195",
//            "194",
//            "193",
//            "192",
//            "191",
//            "190",
//            "189",
//            "188",
//            "187",
//            "186",
//            "185",
//            "184",
//            "183",
//            "182",
//            "181",
//            "180",
//            "179",
//            "178",
//            "177",
//            "176",
//            "175",
//            "174",
//            "173",
//            "172",
//            "171",
//            "170",
//            "169",
//            "168",
//            "167",
//            "166",
//            "165",
//            "164",
//            "163",
//            "162",
//            "161",
//            "160",
//            "159",
//            "158",
//            "157",
//            "156",
//            "155",
//            "154",
//            "153",
//            "152",
//            "151",
//            "150",
//            "149",
//            "148",
//            "147",
//            "146",
//            "145",
//            "144",
//            "143",
//            "142",
//            "141",
//            "140",
//            "139",
//            "138",
//            "137",
//            "136",
//            "135",
//            "134",
//            "133",
//            "132",
//            "131",
//            "130",
//            "129",
//            "128",
//            "127",
//            "126",
//            "125",
//            "124",
//            "123",
//            "122",
//            "121",
//            "120",
//            "119",
//            "118",
//            "117",
//            "116",
//            "115",
//            "114",
//            "113",
//            "112",
//            "111",
//            "110",
//            "109",
//            "108",
//            "107",
//            "106",
//            "105",
//            "104",
//            "103",
//            "102",
//            "101",
//            "100",
//            "99",
//            "98",
//            "97",
//            "96",
//            "95",
//            "94",
//            "93",
//            "92",
//            "91",
//            "90",
//            "89",
//            "88",
//            "87",
//            "86",
//            "85",
//            "84",
//            "83",
//            "82",
//            "81",
//            "80",
//            "79",
//            "78",
//            "77",
//            "76",
//            "75",
//            "74",
//            "73",
//            "72",
//            "71",
//            "70",
//            "69",
//            "68",
//            "67",
//            "66",
//            "65",
//            "64",
//            "63",
//            "62",
//            "61",
//            "60",
//            "59",
//            "58",
//            "57",
//            "56",
//            "55",
//            "54",
//            "53",
//            "52",
//            "51",
//            "50",
//            "49",
//            "48",
//            "47",
//            "46",
//            "45",
//            "44",
//            "43",
//            "42",
//            "41",
//            "40",
//            "39",
//            "38",
//            "37",
//            "36",
//            "35",
//            "34",
//            "33",
//            "32",
//            "31",
//            "30",
//            "29",
//            "28",
//            "27",
//            "26",
//            "25",
//            "24",
//            "23",
//            "22",
//            "21",
//            "20",
//            "19",
//            "18",
//            "17",
//            "16",
//            "15",
//            "14",
//            "13",
//            "12",
//            "11",
//            "10",
//            "9",
//            "8",
//            "7",
//            "6",
//            "5",
//            "4",
//            "3",
//            "2",
//            "1",
           })
    private int enumeration;
    @Param({"true", "false"})
    private boolean environmentFiltering;
    private StructuralMotif queryMotif;

    public static void main(String[] args) throws RunnerException, IOException {
        Options opt = new OptionsBuilder()
                .include(Fit3DBenchmarkEnvironmentFiltering.class.getSimpleName())
                .warmupIterations(5)
                .measurementIterations(5)
                .forks(1)
                .mode(Mode.AverageTime)
//                .timeout(TimeValue.hours(24))
                .timeUnit(TimeUnit.MILLISECONDS)
                .build();
        Collection<RunResult> results = new Runner(opt).run();
        StringJoiner stringJoiner = new StringJoiner("\n", "enumeration,environment_filtering,match_count,size,extent,label_count,min,max,mean,stdev,ci95_min,ci95_max\n", "");
        for (RunResult result : results) {
            Statistics statistics = result.getPrimaryResult().getStatistics();
            double[] confidenceInterval = statistics.getConfidenceIntervalAt(0.95);

            // parse the original motif again to determine extent and size
            String enumeration = result.getParams().getParam("enumeration");
            String environmentFiltering = result.getParams().getParam("environmentFiltering");
            Path queryMotifPath = Paths.get(Resources.getResourceAsFileLocation("structural_motifs/enumerated/" + enumeration + ".pdb"));
            StructuralMotif queryMotif = StructuralMotif.fromLeafSubstructures(StructureParser.local()
                                                                                              .path(queryMotifPath)
                                                                                              .parse().getAllLeafSubstructures());
            double extent = Structures.calculateExtent(queryMotif);
            int labelCount = (int) queryMotif.getAllLeafSubstructures().stream()
                                             .map(LeafSubstructure::getFamily)
                                             .distinct()
                                             .count();

            String resultLine = enumeration +
                                "," +
                                environmentFiltering +
                                "," +
                                queryMotif.size() +
                                "," +
                                extent +
                                "," +
                                labelCount +
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
        Files.write(Paths.get("results_environment_filtering.csv"), stringJoiner.toString().getBytes());
    }

    @Setup
    public void setUp() {
        queryMotif = StructuralMotif.fromLeafSubstructures(StructureParser.local()
                                                                          .path(Paths.get(Resources.getResourceAsFileLocation("structural_motifs/enumerated/" + enumeration + ".pdb")))
                                                                          .parse()
                                                                          .getAllLeafSubstructures());
    }

    @Benchmark
    @Timeout(time = 30, timeUnit = TimeUnit.MINUTES)
    public void runWithLocalMmtf() {
        StructureParser.MultiParser multiParser = StructureParser.local()
                                                                 .localPDB(new StructureParser.LocalPDB(BenchmarkConstants.LOCAL_PDB_LOCATION, SourceLocation.OFFLINE_MMTF))
                                                                 .chainList(BenchmarkConstants.CHAIN_LIST_PATH_10, "\t")
                                                                 .setOptions(BenchmarkConstants.STRUCTURE_PARSER_OPTIONS);
        Fit3DBuilder.ParameterStep parameterStep = Fit3DBuilder.create()
                                                               .query(queryMotif)
                                                               .targets(multiParser)
                                                               .maximalParallelism()
                                                               .atomFilter(StructuralEntityFilter.AtomFilter.isArbitrary());
        Fit3D run;
        if (environmentFiltering) {
            run = parameterStep.filterEnvironments(8.0).run();
        } else {
            run = parameterStep.run();
        }
    }
}
