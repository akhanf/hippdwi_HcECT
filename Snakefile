subjects = ['103']

sessions = ['001']
hemis = ['R','Lnoflip']
metrics = ['MD','FA']
configfile: 'config.yml'

import pandas as pd
df = pd.read_table('participants.tsv')
subjects = df.participant_id.to_list() 
subjects = [ s.strip('sub-') for s in subjects ]


rule all:
    input:
#        expand('results/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_from-dwi_to-Hipp{hemi}.nii.gz',subject=subjects, session=sessions, hemi=hemis),
#        expand('results/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_dwi_proc-FSL_space-Hipp{hemi}_{metric}.nii.gz',subject=subjects, session=sessions, hemi=hemis,metric=metrics)
        expand('results/table_ses-{session}_desc-subfields_{metric}.csv',session=sessions,metric=metrics)

rule rigid_reg_with_init:
    input:
        ref = os.path.join(config['hippunfold_dir'],'sub-{subject}/ses-{session}/hemi-{hemi}/img.nii.gz'),
        flo = os.path.join(config['prepdwi_dir'],'sub-{subject}/ses-{session}/dwi/sub-{subject}_ses-{session}_dwi_proc-FSL_S0.nii.gz'),
        init_xfm = os.path.join(config['hippunfold_dir'],'sub-{subject}/ses-{session}/sub2coronalOblique.txt'),
    params:
        out_prefix = 'results/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_from-dwi_to-Hipp{hemi}_',
#        multires = '--convergence [100x50,1e-6,10]  --shrink-factors 4x2  --smoothing-sigmas 4x1vox'
        multires = '--convergence [40x20x10,1e-6,10]  --shrink-factors 8x4x2  --smoothing-sigmas 8x4x2vox'
    output:
        warped = 'results/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_from-dwi_to-Hipp{hemi}.nii.gz',
        xfm = 'results/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_from-dwi_to-Hipp{hemi}_0GenericAffine.mat',
    container: config['singularity']['ants']
    threads: 8
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsRegistration -d 3 --interpolation Linear {params.multires} --metric CC[{input.ref},{input.flo},1,3] --transform Rigid[0.1] --output [{params.out_prefix},{output.warped}] --initial-moving-transform [{input.init_xfm}]  -v '



rule apply_xfm_to_dti_metrics:
    input:
        ref = os.path.join(config['hippunfold_dir'],'sub-{subject}/ses-{session}/hemi-{hemi}/img.nii.gz'),
        flo = os.path.join(config['prepdwi_dir'],'sub-{subject}/ses-{session}/dwi/sub-{subject}_ses-{session}_dwi_proc-FSL_{metric}.nii.gz'),
        xfm = 'results/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_from-dwi_to-Hipp{hemi}_0GenericAffine.mat',
    output: 
        warped = 'results/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_dwi_proc-FSL_space-Hipp{hemi}_{metric}.nii.gz'
    threads: 8
    container: config['singularity']['ants']
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation Linear -i {input.flo} -o {output.warped} -r {input.ref} -t {input.xfm}'

rule gen_mean_value_txt:
    input:
        dti = 'results/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_dwi_proc-FSL_space-Hipp{hemi}_{metric}.nii.gz',
        subfields = os.path.join(config['hippunfold_dir'],'sub-{subject}/ses-{session}/hemi-{hemi}/subfields-BigBrain.nii.gz')
    params:
        hemi_letter = lambda wildcards: wildcards.hemi[0]
    output:
        csv = 'results/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_hemi-{hemi}_desc-subfields_{metric}.csv'
    container: config['singularity']['fsl']
    shell:
        "fslstats -K {input.subfields} {input.dti} -m  > {output.csv} "

rule concat_lr_csv:
    input:
        left = 'results/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_hemi-Lnoflip_desc-subfields_{metric}.csv',
        right = 'results/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_hemi-R_desc-subfields_{metric}.csv'
    output: 'results/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_desc-subfields_{metric}.csv'
    shell:
        "left=`cat {input.left}` && right=`cat {input.right}` && echo sub-{wildcards.subject} $left $right | sed 's/ /,/g'> {output}"

rule create_group_csv:
    input: 
        expand('results/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_desc-subfields_{metric}.csv',subject=subjects,allow_missing=True)
    params:
        header = 'subject,L_SUB,L_CA1,L_CA2,L_CA3,L_CA4,L_DG,R_SUB,R_CA1,R_CA2,R_CA3,R_CA3,R_DG'
    output: 'results/table_ses-{session}_desc-subfields_{metric}.csv'
    shell:
        'echo {params.header} > {output} && cat {input} >> {output}'

        
