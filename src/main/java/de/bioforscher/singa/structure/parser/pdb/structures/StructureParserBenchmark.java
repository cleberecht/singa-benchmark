package de.bioforscher.singa.structure.parser.pdb.structures;

import de.bioforscher.singa.structure.BenchmarkConstants;
import de.bioforscher.singa.structure.parser.pdb.structures.StructureParser.LocalPDB;
import de.bioforscher.singa.structure.parser.pdb.structures.StructureParser.MultiParser;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;
import org.openjdk.jmh.runner.options.TimeValue;

import java.util.concurrent.TimeUnit;

/**
 * @author fk
 */
@State(Scope.Benchmark)
public class StructureParserBenchmark {

    private StructureParserOptions structureParserOptions;

    public static void main(String[] args) throws RunnerException {
        Options opt = new OptionsBuilder()
                .include(StructureParserBenchmark.class.getSimpleName())
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
        structureParserOptions = StructureParserOptions.withSettings(
                StructureParserOptions.Setting.OMIT_EDGES,
                StructureParserOptions.Setting.OMIT_HYDROGENS,
                StructureParserOptions.Setting.OMIT_LIGAND_INFORMATION,
                StructureParserOptions.Setting.GET_IDENTIFIER_FROM_FILENAME);
    }


    @Benchmark
    public void batchLocalPdb() {
        MultiParser multiParser = StructureParser.local()
                                                 .localPDB(new LocalPDB(BenchmarkConstants.LOCAL_PDB_LOCATION, SourceLocation.OFFLINE_PDB))
                                                 .chainList(BenchmarkConstants.CHAIN_LIST_PATH_500, "\t")
                                                 .setOptions(structureParserOptions);
        while (multiParser.hasNext()) {
            try {
                multiParser.next().getAllAtoms();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    @Benchmark
    public void batchLocalMmtf() {
        MultiParser multiParser = StructureParser.local()
                                                 .localPDB(new LocalPDB(BenchmarkConstants.LOCAL_PDB_LOCATION, SourceLocation.OFFLINE_MMTF))
                                                 .chainList(BenchmarkConstants.CHAIN_LIST_PATH_500, "\t")
                                                 .setOptions(structureParserOptions);
        while (multiParser.hasNext()) {
            try {
                multiParser.next().getAllAtoms();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
}
