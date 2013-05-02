package org.rostlab.master.practical2013.alignment;

import java.util.Arrays;

import org.rostlab.master.practical2013.RunConfiguration;

/**
 * This class provides user interface for querying alignment results
 * @author Shen Wei
 *
 */
public class ResultViewer {
	
	public static void main(String[] args) {
		
		try {
			RunConfiguration rc = RunConfiguration.load(args[0]);
			String pathResult = args[1];
			String command = args[2];
			
			if (command.equals("extract")) {
				cmdExtract(
						rc, 
						pathResult, 
						Arrays.copyOfRange(args, 3, args.length-1,args.getClass()));
			}
			
		} catch (ArrayIndexOutOfBoundsException e) {
			System.err.println("View alignment results\n" +
					"Usage: result-viewer <run_configuration.conf> <result_file> <command> [parameters]");
			System.err.println("Commands:\n" +
					"\textract <seq_id> [start_pos_inclusive [end_pos_inclusive]]\n" +
					"\tbatch <command_table_file>\n");
		}
	}

	private static void cmdExtract(RunConfiguration rc, String pathResult,
			String[] parameters) {
		// TODO Auto-generated method stub
		
	}
}
