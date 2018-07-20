import bio.singa.structure.model.oak.StructuralMotif;
import bio.singa.structure.model.oak.Structures;
import bio.singa.structure.parser.pdb.structures.StructureParser;
import bio.singa.structure.parser.pdb.structures.StructureWriter;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author fk
 */
public class DataPreparation {

    public static final Path OUTPUT_PATH = Paths.get("/home/fkaiser/Workspace/IdeaProjects/singa-benchmark/src/main/resources/structural_motifs");
    private static final Logger logger = LoggerFactory.getLogger(DataPreparation.class);
    private static final double MOTIF_EXTENT_CUTOFF = 15.0;
    public static final Predicate<StructuralMotif> MOTIF_FILTER = structuralMotif ->
//            structuralMotif.getAllLeafSubstructures().stream()
//                           .map(LeafSubstructure::getFamily)
//                           .distinct()
//                           .count() == structuralMotif.size() &&
            Structures.calculateExtent(structuralMotif) <= MOTIF_EXTENT_CUTOFF;
    public static final int MINIMAL_MOTIF_SIZE = 2;
    public static final int MAXIMAL_MOTIF_SIZE = 6;

    @Test
    public void groupCsaMotifsBySize() throws IOException {

        Path csaMotifPath = Paths.get("/home/fkaiser/Workspace/CloudStation/PhD/Promotion/datasets/csa/motifs_BLAST_10e80");
        List<Path> motifPaths = Files.walk(csaMotifPath)
                                     .filter(path -> path.toFile().isFile())
                                     .collect(Collectors.toList());

        logger.info("handling {} motifs from CSA", motifPaths.size());

        Map<Integer, List<StructuralMotif>> motifsBySize = new TreeMap<>();
        for (Path motifPath : motifPaths) {
            logger.info("parsing motif {}", motifPath);
            StructuralMotif structuralMotif = StructuralMotif.fromLeafSubstructures(StructureParser.local()
                                                                                                   .path(motifPath)
                                                                                                   .parse().getAllLeafSubstructures());
            if (!MOTIF_FILTER.test(structuralMotif)|| structuralMotif.size() > MAXIMAL_MOTIF_SIZE || structuralMotif.size() < MINIMAL_MOTIF_SIZE) {
                logger.info("skipping motif {} that does not meet filter criteria", motifPath);
                continue;
            }
            int size = structuralMotif.size();
            if (motifsBySize.containsKey(size)) {
                motifsBySize.get(size).add(structuralMotif);
            } else {
                motifsBySize.put(size, Stream.of(structuralMotif).collect(Collectors.toList()));
            }
        }

        for (Map.Entry<Integer, List<StructuralMotif>> entry : motifsBySize.entrySet()) {
            for (StructuralMotif structuralMotif : entry.getValue()) {
                StructureWriter.writeLeafSubstructureContainer(structuralMotif, OUTPUT_PATH.resolve("size").resolve("size_" + entry.getKey()).resolve(structuralMotif.toString() + ".pdb"));
            }
        }
    }

    @Test
    public void groupCsaMotifsByExtent() throws IOException {
        Path csaMotifPath = Paths.get("/home/fkaiser/Workspace/CloudStation/PhD/Promotion/datasets/csa/motifs_BLAST_10e80");
        List<Path> motifPaths = Files.walk(csaMotifPath)
                                     .filter(path -> path.toFile().isFile())
                                     .collect(Collectors.toList());

        logger.info("handling {} motifs from CSA", motifPaths.size());

        Map<Integer, List<StructuralMotif>> motifsByExtent = new TreeMap<>();
        for (Path motifPath : motifPaths) {
            logger.info("parsing motif {}", motifPath);
            StructuralMotif structuralMotif = StructuralMotif.fromLeafSubstructures(StructureParser.local()
                                                                                                   .path(motifPath)
                                                                                                   .parse().getAllLeafSubstructures());

            if (!MOTIF_FILTER.test(structuralMotif)|| structuralMotif.size() > MAXIMAL_MOTIF_SIZE || structuralMotif.size() < MINIMAL_MOTIF_SIZE) {
                logger.info("skipping motif {} that does not meet filter criteria", motifPath);
                continue;
            }

            double extent = Structures.calculateExtent(structuralMotif);
            int bin = (int) Math.floor(extent);
            if (motifsByExtent.containsKey(bin)) {
                motifsByExtent.get(bin).add(structuralMotif);
            } else {
                motifsByExtent.put(bin, Stream.of(structuralMotif).collect(Collectors.toList()));
            }
        }

        for (Map.Entry<Integer, List<StructuralMotif>> entry : motifsByExtent.entrySet()) {
            for (StructuralMotif structuralMotif : entry.getValue()) {
                StructureWriter.writeLeafSubstructureContainer(structuralMotif, OUTPUT_PATH.resolve("extent").resolve("extent_" + entry.getKey()).resolve(structuralMotif.toString() + ".pdb"));
            }
        }
    }

    @Test
    public void enumerateCsaMotifs() throws IOException {

        Path csaMotifPath = Paths.get("/home/fkaiser/Workspace/CloudStation/PhD/Promotion/datasets/csa/motifs_BLAST_10e80");

        List<Path> motifPaths = Files.walk(csaMotifPath)
                                     .filter(path -> path.toFile().isFile())
                                     .collect(Collectors.toList());

        logger.info("handling {} motifs from CSA", motifPaths.size());

        List<StructuralMotif> structuralMotifs = new ArrayList<>();
        for (Path motifPath : motifPaths) {
            logger.info("parsing motif {}", motifPath);
            StructuralMotif structuralMotif = StructuralMotif.fromLeafSubstructures(StructureParser.local()
                                                                                                   .path(motifPath)
                                                                                                   .parse().getAllLeafSubstructures());

            if (!MOTIF_FILTER.test(structuralMotif) || structuralMotif.size() > MAXIMAL_MOTIF_SIZE || structuralMotif.size() < MINIMAL_MOTIF_SIZE) {
                logger.info("skipping motif {} that does not meet filter criteria", motifPath);
                continue;
            }
            structuralMotifs.add(structuralMotif);
        }
        structuralMotifs.sort(Comparator.comparing(StructuralMotif::size)
                                        .thenComparing(Structures::calculateExtent));
        for (int i = 0; i < structuralMotifs.size(); i++) {
            StructuralMotif structuralMotif = structuralMotifs.get(i);
            StructureWriter.writeLeafSubstructureContainer(structuralMotif, OUTPUT_PATH.resolve("enumerated").resolve((i + 1) + ".pdb"));
        }
    }
}
