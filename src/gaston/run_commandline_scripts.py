import os
import subprocess

def train_NN_parallel(path_to_coords, path_to_glmpca, hidden_spatial, hidden_expression, 
                      output_dir, conda_environment, path_to_conda, pos_encoding=True,
                      epochs=10000, checkpoint=500, optimizer='adam', num_seeds=30, parallel_processes=4):

    hidden_spatial = ' '.join(map(str, hidden_spatial))
    hidden_expression = ' '.join(map(str, hidden_expression))
    
    commands = []
    for seed in range(num_seeds):
        os.makedirs(f"{output_dir}/seed{seed}", exist_ok=True)
        cmd = f"gaston -i {path_to_coords} -o {path_to_glmpca} "
        cmd += f"--epochs {epochs} -d {output_dir} "
        cmd += f"--hidden_spatial {hidden_spatial} --hidden_expression {hidden_expression} "
        cmd += f"--optimizer {optimizer} --seed {seed} -c {checkpoint}"
        if pos_encoding:
            cmd += f" --positional_encoding"
        commands.append(cmd)

    # Create the bash script
    script_path = f"{output_dir}/run_all_seeds.sh"
    with open(script_path, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("echo 'Starting training for all seeds in parallel...'\n\n")
        
        # Add conda environment activation to the script
        f.write(f"source {path_to_conda} {conda_environment} \n \n")

        # Add the commands to the script and run them in parallel
        commands_str = '\n'.join(commands)
        f.write(f"echo '{commands_str}' | xargs -P {parallel_processes} -I CMD bash -c CMD\n")

    # Make the script executable
    os.chmod(script_path, 0o755)

    # Run the created script
    subprocess.run([script_path])

    os.remove(script_path)

    print(f"Bash script executed: {script_path}")