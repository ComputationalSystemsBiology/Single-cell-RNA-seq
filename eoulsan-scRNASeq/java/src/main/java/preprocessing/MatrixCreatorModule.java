package preprocessing;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import fr.ens.biologie.genomique.eoulsan.Globals;
import fr.ens.biologie.genomique.eoulsan.annotations.LocalOnly;
import fr.ens.biologie.genomique.eoulsan.core.InputPorts;
import fr.ens.biologie.genomique.eoulsan.core.InputPortsBuilder;
import fr.ens.biologie.genomique.eoulsan.core.OutputPorts;
import fr.ens.biologie.genomique.eoulsan.core.TaskContext;
import fr.ens.biologie.genomique.eoulsan.core.TaskResult;
import fr.ens.biologie.genomique.eoulsan.core.TaskStatus;
import fr.ens.biologie.genomique.eoulsan.core.Version;
import fr.ens.biologie.genomique.eoulsan.data.DataFormat;
import fr.ens.biologie.genomique.eoulsan.data.DataFormatRegistry;
import fr.ens.biologie.genomique.eoulsan.modules.AbstractModule;
import fr.ens.biologie.genomique.eoulsan.data.Data;
import fr.ens.biologie.genomique.eoulsan.data.DataFormats;


import static fr.ens.biologie.genomique.eoulsan.EoulsanLogger.getLogger;
import static fr.ens.biologie.genomique.eoulsan.core.OutputPortsBuilder.singleOutputPort;

/**
 * This class define a module that merges expression files.
 *
 * @author Geoffray Brelurut
 * @since 2017
 */

@LocalOnly
public class MatrixCreatorModule extends AbstractModule {
    /**
     * Module Name
     */
    private static final String MODULE_NAME = "matrixcreator";

    // DataFormat
    private static DataFormat EXPRESSION_MATRIX_TSV =
            DataFormatRegistry.getInstance().getDataFormatFromName("expression_matrix_tsv");

    //
    // Module Methods
    //

    @Override
    public String getName() {
        return MODULE_NAME;
    }

    @Override
    public Version getVersion() {
        return Globals.APP_VERSION;
    }

    @Override
    public String getDescription() {
        return "This class merges the files from expression step";
    }

    @Override
    public InputPorts getInputPorts() {
        return new InputPortsBuilder()
                .addPort(InputPortsBuilder.DEFAULT_SINGLE_INPUT_PORT_NAME, true, DataFormats.EXPRESSION_RESULTS_TSV)
                .create();
    }

    @Override
    public OutputPorts getOutputPorts() {
        return singleOutputPort("matrix", EXPRESSION_MATRIX_TSV);
    }

    @Override
    public TaskResult execute(TaskContext context, TaskStatus status) {

        // Get and Check Input
        final Data inData = context.getInputData(DataFormats.EXPRESSION_RESULTS_TSV);

        if (!inData.isList() || inData.getListElements().size() < 2) {
            getLogger().severe("Only one expression file found");
            return status.createTaskResult();
        }

        // Get output and input files
        File output = context.getOutputData(EXPRESSION_MATRIX_TSV, inData).getDataFile().toFile();
        //String names = getSamplesName(context);
        List<File> files = new ArrayList<>();
        for (Data data : inData.getListElements()) {
            File file = data.getDataFile().toFile();
            files.add(file);
        }

        // Merge files
        try {
            merge(files, output);
            return status.createTaskResult();
        } catch (NullPointerException | IOException e) {
            return status.createTaskResult(e, "Error with file :" + e.getMessage());
        }
    }

    /**
     * Merge two columns files into one matrix
     *
     * @param files      list of files to merge
     * @param mergedFile file for writing merging
     *                   //@param headers headers of the matrix
     * @throws NullPointerException
     * @throws IOException
     */

    private static void merge(List<File> files, File mergedFile)
            throws NullPointerException, IOException {

			/* Test parameters */
        if (files == null) {
            throw new NullPointerException("files argument cannot be null");
        }

        if (mergedFile == null) {
            throw new NullPointerException("mergedFile argument cannot be null");
        }

			/* Set structures and extract first file */
        File init = files.get(0);

        String headers = init.getName().replace(".tsv", "").split("_")[3];
        List<String> genes = new ArrayList<>();
        Map<String, String> lines = new HashMap<>();

			/* Treat first file initialising structures */
        try (BufferedReader in = new BufferedReader(new FileReader(init))) {
            String aLine;
            String[] parts = new String[2];
            int i = 0;
            while ((aLine = in.readLine()) != null) {
                if (i > 0) {
                    parts = aLine.split("\t");
                    genes.add(parts[0]);
                    lines.put(parts[0], parts[1]);
                }
                i++;
            }
        }

			/* Remove file treated from files list */
        files.remove(0);

			/* Treat other files */
        for (File file : files) {

            headers = headers + "\t" + file.getName().replace(".tsv", "").split("_")[3];
            try (BufferedReader in = new BufferedReader(new FileReader(file))) {
                String aLine;
                String[] parts = new String[2];
                while ((aLine = in.readLine()) != null) {
                    parts = aLine.split("\t");
                    String value = lines.get(parts[0]);
                    value = value + "\t" + parts[1];
                    lines.put(parts[0], value);
                }
            }
        }

			/* Write merged file, excluding non detected genes */
        try (BufferedWriter out = new BufferedWriter(new FileWriter(mergedFile,
                true))) {
            out.write(headers);
            out.newLine();
            for (String key : genes) {
                int sum = 0;
                String[] values = lines.get(key).split("\t");
                for (String value : values) {
                    sum = sum + Integer.parseInt(value);
                }
                if (sum > 0) {
                    out.write(key + "\t" + lines.get(key));
                    out.newLine();
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Extract samples' name from context, returning them in tabulated string
     * @param context TaskContext of the module
     * @return String gathering samples' name (for tabulated file headers)
     */
    /*private static String getSamplesName(TaskContext context) {
		Design design = context.getWorkflow().getDesign();
		List<Sample> samples = design.getSamples();
		StringBuilder names = new StringBuilder();
		for(Sample sample : samples) {
			names.append("\t").append(sample.getName());
		}
		return names.toString();
	}*/

}
