package fr.ens.biologie.genomique.eoulsan.modules.preprocessing;

import java.io.*;
import java.nio.file.Files;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import fr.ens.biologie.genomique.eoulsan.annotations.LocalOnly;
import fr.ens.biologie.genomique.eoulsan.data.DataFormat;
import fr.ens.biologie.genomique.eoulsan.data.DataFormatRegistry;
import fr.ens.biologie.genomique.eoulsan.design.io.DesignFormatFinderInputStream;
import fr.ens.biologie.genomique.eoulsan.modules.AbstractModule;
import fr.ens.biologie.genomique.eoulsan.Globals;
import fr.ens.biologie.genomique.eoulsan.core.OutputPorts;
import fr.ens.biologie.genomique.eoulsan.core.TaskContext;
import fr.ens.biologie.genomique.eoulsan.core.TaskResult;
import fr.ens.biologie.genomique.eoulsan.core.TaskStatus;
import fr.ens.biologie.genomique.eoulsan.core.Version;
import fr.ens.biologie.genomique.eoulsan.data.DataFile;

import static fr.ens.biologie.genomique.eoulsan.core.OutputPortsBuilder.singleOutputPort;


/**
 * This class define a module that extract columns section from design file
 *
 * @author Geoffray Brelurut
 * @since 2017
 */

@LocalOnly
public class CellDataExtractorModule extends AbstractModule {

    /**
     * Module name
     */
    private static final String MODULE_NAME = "reducer";

    // DataFormat
    private static DataFormat INI_CELLS_METADATA_TSV =
            DataFormatRegistry.getInstance().getDataFormatFromName("initial_cells_metadata_tsv");

    //
    // Module methods
    //

    @Override
    public String getName() {
        return MODULE_NAME;
    }

    @Override
    public String getDescription() {
        return "This module extracts columns from design file to further use in scRNA-Seq analysis";
    }

    @Override
    public Version getVersion() {
        return Globals.APP_VERSION;
    }


    @Override
    public OutputPorts getOutputPorts() {
        return singleOutputPort("cellsoutput", INI_CELLS_METADATA_TSV);
    }

    @Override
    public TaskResult execute(TaskContext context, TaskStatus status) {

        // Get input
        final DataFile design = context.getDesignFile();
        final File input = design.toFile();


        // Get output
        final File output = context.getOutputData(INI_CELLS_METADATA_TSV, design.getName().replace(".txt", "")).getDataFile().toFile();

        // Process file
        try {
            // Test design format version
            DesignFormatFinderInputStream formatFinder = new DesignFormatFinderInputStream(new FileInputStream(input));
            int version = formatFinder.getDesignFormatVersion();

            // If design v2 extract table
            if ( version == 2) {
                reduceDesign(input, output);
            } else { // Else make a copy
                Files.copy(input.toPath(), output.toPath());
            }
            return status.createTaskResult();
        } catch (final IOException e) {
            return status.createTaskResult(e, "Error with file :" + e.getMessage());
        }
    }

    /**
     * Parse design file extracting columns section
     *
     * @param file   the design file
     * @param output the output file
     * @throws IOException if an error occurs while parsing or writing file
     */
    protected static void reduceDesign(File file, File output) throws IOException {
        /* Initialise variables */
        String aLine;
        boolean table = false;
        String pattern = "(Columns)";
        Pattern seek = Pattern.compile(pattern);

		/* Read file, copy lines */
        try (BufferedReader in = new BufferedReader(new FileReader(file))) {
            while ((aLine = in.readLine()) != null) {
                if (table) {
                    try (BufferedWriter out = new BufferedWriter(new FileWriter(output, true))) {
                        out.write(aLine);
                        out.newLine();
                    }
                } else {
                    Matcher m = seek.matcher(aLine);
                    if (m.find()) {
                        table = true;
                    }
                }
            }
        }
    }
}
