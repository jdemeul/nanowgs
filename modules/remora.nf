
/* 
* Remora modified basecall training – TO BE IMPLEMENTED
*/
process remora {
    label 'remora'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path sorted_bam

    output:
    path "snv_indel_medaka", emit: medaka_variant

    script:
    """
        # PREPARE INPUT DATA USING MEGA AND TAIYAKI

        singularity run --nv -B /staging/leuven/stg_00002/lcb/jdemeul/ \
            /staging/leuven/stg_00002/lcb/jdemeul/software/singularity_images/zeunas-guppy-5.0.16-megalodon_v2.4.1-rerio-taiyaki_v5.3.0.img \
            megalodon \
            ./Unmod/ \
            --reference \
            /staging/leuven/stg_00002/lcb/jdemeul/reference/chm13_v1.1_chrY_KI270740_EBV/indexes/minimap2-ont/genome.mmi \
            --output-directory mega_res_pcr \
            --outputs mappings signal_mappings \
            --num-reads 10000 \
            --guppy-config dna_r9.4.1_450bps_sup.cfg \
            --processes 20 \
            --guppy-server-path /opt/ont-guppy/bin/guppy_basecall_server \
            --devices "cuda:all" \
            --overwrite

        singularity run --nv -B /staging/leuven/stg_00002/lcb/jdemeul/ \
            /staging/leuven/stg_00002/lcb/jdemeul/software/singularity_images/zeunas-guppy-5.0.16-megalodon_v2.4.1-rerio-taiyaki_v5.3.0.img \
            megalodon \
            ./Hia5/ \
            --reference /staging/leuven/stg_00002/lcb/jdemeul/reference/chm13_v1.1_chrY_KI270740_EBV/indexes/minimap2-ont/genome.mmi \
            --output-directory mega_res_hia5 \
            --ref-mods-all-motifs Y 6mA A 0 \
            --outputs mappings signal_mappings \
            --num-reads 10000 \
            --guppy-config dna_r9.4.1_450bps_sup.cfg \
            --processes 20 \
            --guppy-server-path /opt/ont-guppy/bin/guppy_basecall_server \
            --devices "cuda:all" \
            --overwrite


        # MERGE TRAINING DATA SETS USING TAIYAKI

        singularity run --nv -B /staging/leuven/stg_00002/lcb/jdemeul/ \
            /staging/leuven/stg_00002/lcb/jdemeul/software/singularity_images/zeunas-guppy-5.0.16-megalodon_v2.4.1-rerio-taiyaki_v5.3.0.img \
            python3 /taiyaki/misc/merge_mappedsignalfiles.py mapped_signal_train_data.hdf5 \
            --input mega_res_pcr/signal_mappings.hdf5 None \
            --input mega_res_hia5/signal_mappings.hdf5 None \
            --allow_mod_merge \
            --batch_format
        # Merged alphabet contains: canonical alphabet ACGT with modified base(s) Y=6mA (alt to A)
        # Writing reads to mapped_signal_train_data.hdf5
        # Copied 7077 reads from mega_res_pcr/signal_mappings.hdf5.
        # Copied 4296 reads from mega_res_hia5/signal_mappings.hdf5.
        # Copied 11373 reads in total.


        # REMORA DATA PREP

        singularity run --nv -B /staging/leuven/stg_00002/lcb/jdemeul/ \
            /staging/leuven/stg_00002/lcb/jdemeul/software/singularity_images/zeunas-remora-0.1.1-569ded0.img \
            remora dataset prepare mapped_signal_train_data.hdf5 \
            --output-remora-training-file remora_train_chunks.npz \
            --motif A 0 \
            --mod-bases Y \
            --chunk-context 50 50 \
            --kmer-context-bases 6 6 \
            --max-chunks-per-read 20 \
            --log-filename log.txt
        # [21:51:38] Opening mapped signal files
        # [21:51:38] Allocating memory for output tensors
        # [21:51:38] Processing reads
        # 100%|████████████████████████████████████████████████████████████████████████████████████████| 11373/11373 [00:27<00:00, 409.38reads/s]
        # [21:52:06] Processing encountered 0 invalid reads from Taiyaki and 1116 short chunks which were discarded.
        # [21:52:06] Extracted 226279 chunks from 11373 reads.
        # [21:52:06] Label distribution: Counter({0: 140806, 1: 85473})

        # REMORA TRAINING

        singularity run --nv -B /staging/leuven/stg_00002/lcb/jdemeul/ \
            /staging/leuven/stg_00002/lcb/jdemeul/software/singularity_images/zeunas-remora-0.1.1-569ded0.img \
            remora model train remora_train_chunks.npz \
            --model remora/models/Conv_w_ref.py \
            --device 0 \
            --output-path remora_train_results

        # [22:41:00] Seed selected is 844325128
        # [22:41:03] Loading dataset from Remora file
        # [22:41:03] Loaded data info from file:
        #           base_pred : False
        #           mod_bases : Y
        #  kmer_context_bases : (6, 6)
        #       chunk_context : (50, 50)
        #               motif : ('A', 0)

        # [22:41:03] Loading model
        # [22:41:03] Model structure:
        # network(
        #   (sig_conv1): Conv1d(1, 4, kernel_size=(11,), stride=(1,))
        #   (sig_conv2): Conv1d(4, 16, kernel_size=(11,), stride=(1,))
        #   (sig_conv3): Conv1d(16, 64, kernel_size=(9,), stride=(3,))
        #   (seq_conv1): Conv1d(52, 16, kernel_size=(11,), stride=(1,))
        #   (seq_conv2): Conv1d(16, 32, kernel_size=(11,), stride=(1,))
        #   (seq_conv3): Conv1d(32, 64, kernel_size=(9,), stride=(3,))
        #   (merge_conv1): Conv1d(128, 64, kernel_size=(5,), stride=(1,))
        #   (merge_conv2): Conv1d(64, 64, kernel_size=(5,), stride=(1,))
        #   (merge_conv3): Conv1d(64, 64, kernel_size=(3,), stride=(2,))
        #   (merge_conv4): Conv1d(64, 64, kernel_size=(3,), stride=(2,))
        #   (fc): Linear(in_features=192, out_features=2, bias=True)
        #   (sig_bn1): BatchNorm1d(4, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
        #   (sig_bn2): BatchNorm1d(16, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
        #   (sig_bn3): BatchNorm1d(64, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
        #   (seq_bn1): BatchNorm1d(16, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
        #   (seq_bn2): BatchNorm1d(32, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
        #   (seq_bn3): BatchNorm1d(64, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
        #   (merge_bn1): BatchNorm1d(64, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
        #   (merge_bn2): BatchNorm1d(64, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
        #   (merge_bn3): BatchNorm1d(64, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
        #   (merge_bn4): BatchNorm1d(64, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
        # )
        # [22:41:03]  Params (k) 130.94 | MACs (M) 3303.18
        # [22:41:03] Preparing training settings
        # [22:41:10] Label distribution: Counter({0: 140806, 1: 85473})
        # [22:41:10] Running initial validation
        # [22:41:10] Start training
        # Epochs:  22%|███▌            | 11/50 [03:54<13:50, 21.30s/it, acc_train=0.8086, acc_val=0.7533, loss_train=0.414223, loss_val=0.504485]

    """

}

