#!/usr/bin/env nextflow
import com.google.common.hash.Hashing
import com.google.common.base.Charsets
import groovy.json.JsonOutput

params.help = false
params.input = false

if(params.help) {
    usage = file("$baseDir/USAGE")

    cpu_count = Runtime.runtime.availableProcessors()
    bindings = ["use_provided_centroids":"$params.use_provided_centroids",
                "nb_points":"$params.nb_points",
                "min_lesion_vol":"$params.min_lesion_vol",
                "min_streamline_count":"$params.min_streamline_count",
                "min_voxel_count":"$params.min_voxel_count",
                "mean_std_density_weighting":"$params.mean_std_density_weighting",
                "mean_std_per_point_density_weighting":"$params.mean_std_per_point_density_weighting",
                "endpoints_metric_stats_normalize_weights":"$params.endpoints_metric_stats_normalize_weights",
                "skip_projection_endpoints_metrics":"$params.skip_projection_endpoints_metrics",
                "cpu_count":"$cpu_count",
                "bundle_suffix_to_remove":"$params.bundle_suffix_to_remove",
                "compute_fixel_mrds_metrics":"$params.compute_fixel_mrds_metrics"
        ]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info "Scilpy dMRI tractometry pipeline"
log.info "================================"
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

Channel
    .fromFilePairs("$params.input/**/bundles/*.trk",
                   size: -1) { it.parent.parent.name }
    .into{bundles_for_rm_invalid; in_bundles_check}

Channel
    .fromFilePairs("$params.input/**/metrics/*.nii.gz",
                   size: -1) { it.parent.parent.name }
    .set{in_metrics}

Channel
    .fromFilePairs("$params.input/**/centroids/*.trk",
        size: -1) { it.parent.parent.name }
    .into{centroids_for_resample; in_centroids_check}

Channel
    .fromFilePairs("$params.input/**/*lesion_mask.nii.gz",
        size: -1) { it.parent.name }
    .set{lesion_for_lesion_load}

Channel
    .fromFilePairs("$params.input/**/*fodf.nii.gz",
        size: -1) { it.parent.name }
    .set{fodf_for_fixel_afd}

Channel
    .fromFilePairs("$params.input/**/*pdds.nii.gz",
        size: -1) { it.parent.name }
    .set{pdds_for_fixel_mrds}

Channel
    .fromFilePairs("$params.input/**/fixel_metrics/*fixel_fa.nii.gz",
        size: -1) { it.parent.parent.name }
    .set{fa_for_fixel_mrds}

Channel
    .fromFilePairs("$params.input/**/fixel_metrics/*fixel_md.nii.gz",
        size: -1) { it.parent.parent.name }
    .set{md_for_fixel_mrds}

Channel
    .fromFilePairs("$params.input/**/fixel_metrics/*fixel_rd.nii.gz",
        size: -1) { it.parent.parent.name }
    .set{rd_for_fixel_mrds}

Channel
    .fromFilePairs("$params.input/**/fixel_metrics/*fixel_ad.nii.gz",
        size: -1) { it.parent.parent.name }
    .set{ad_for_fixel_mrds}

in_metrics
    .set{metrics_for_rename}


in_bundles_check.map{it[1]}.flatten().count().set{number_bundles_for_compare}
in_centroids_check.map{it[1]}.flatten().count().set{number_centroids_for_compare}

if (params.use_provided_centroids){
number_centroids_for_compare
    .concat(number_bundles_for_compare)
    .toList()
    .subscribe{a, b -> if (a < b)
    error "Error ~ You ask the pipeline to use provided centroids but there are less centroids than bundles.\nPlease provide at least a centroid per bundle."}
}

process Rename_Metrics {
    input:
    set sid, file(metrics) from metrics_for_rename

    output:
    set sid, "*_metric.nii.gz" into metrics_for_mean_std,
        metrics_for_endpoints_metrics, metrics_for_endpoints_roi_stats,
        metrics_for_volume, metrics_for_mean_std_per_point

    script:
    """
    for metric in *.nii.gz; do
        mv \$metric \$(basename \${metric/${sid}__/} .nii.gz)_metric.nii.gz
    done
    """
}

process Remove_Invalid_Streamlines {
    input:
    set sid, file(bundles) from bundles_for_rm_invalid

    output:
    set sid, "${sid}__*_ic.trk" into bundles_for_label_and_distance_map, bundles_for_centroids, bundles_for_fixel_afd, bundles_for_fixel_mrds

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    for bundle in $bundles_list;
      do if [[ \$bundle == *"__"* ]]; then
          pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
          bname=\${bundle:\$pos}
          bname=\$(basename \$bname .trk)
      else
          bname=\$(basename \$bundle .trk)
      fi
      bname=\${bname/$params.bundle_suffix_to_remove/}

      scil_remove_invalid_streamlines.py \$bundle ${sid}__\${bname}_ic.trk --remove_single_point --remove_overlapping_points --cut_invalid --no_empty
    done
    """
}

bundles_for_fixel_afd
    .join(fodf_for_fixel_afd)
    .set{bundle_fodf_for_fixel_afd}

process Fixel_AFD {
    input:
    set sid, file(bundles), file(fodf) from bundle_fodf_for_fixel_afd

    output:
    set sid, "*_afd_metric.nii.gz" into fixel_afd_for_mean_std,
        fixel_afd_for_endpoints_metrics, fixel_afd_for_endpoints_roi_stats,
        fixel_afd_for_mean_std_per_point

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    for bundle in $bundles_list;
        do if [[ \$bundle == *"__"* ]]; then
            pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
            bname=\${bundle:\$pos}
            bname=\$(basename \$bname .trk)
        else
            bname=\$(basename \$bundle .trk)
        fi
        bname=\${bname/$params.bundle_suffix_to_remove/}
        scil_compute_mean_fixel_afd_from_bundles.py \$bundle $fodf \${bname}_afd_metric.nii.gz
    done
    """
}

bundles_for_fixel_mrds
    .join(pdds_for_fixel_mrds)
    .join(fa_for_fixel_mrds)
    .join(md_for_fixel_mrds)
    .join(rd_for_fixel_mrds)
    .join(ad_for_fixel_mrds)
    .set{bundle_pdds_for_fixel_mrds}

process Fixel_MRDS {
    input:
        tuple sid, file(bundles), file(pdd), file(fa), file(md), file(rd), file(ad) from bundle_pdds_for_fixel_mrds

    output:
        set sid, "*_fa_metric.nii.gz" into fixel_fa_for_mean_std,
            fixel_fa_for_endpoints_metrics, fixel_fa_for_endpoints_roi_stats,
            fixel_fa_for_mean_std_per_point
        set sid, "*_md_metric.nii.gz" into fixel_md_for_mean_std,
            fixel_md_for_endpoints_metrics, fixel_md_for_endpoints_roi_stats,
            fixel_md_for_mean_std_per_point
        set sid, "*_rd_metric.nii.gz" into fixel_rd_for_mean_std,
            fixel_rd_for_endpoints_metrics, fixel_rd_for_endpoints_roi_stats,
            fixel_rd_for_mean_std_per_point
        set sid, "*_ad_metric.nii.gz" into fixel_ad_for_mean_std,
            fixel_ad_for_endpoints_metrics, fixel_ad_for_endpoints_roi_stats,
            fixel_ad_for_mean_std_per_point

    when:
    params.compute_fixel_mrds_metrics

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    for bundle in $bundles_list;
        do if [[ \$bundle == *"__"* ]]; then
            pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
            bname=\${bundle:\$pos}
            bname=\$(basename \$bname .trk)
        else
            bname=\$(basename \$bundle .trk)
        fi
        scil_compute_mean_fixel_lobe_mrds_metric_from_bundles.py \$bundle \
            $pdd $fa \${bname}_fa_metric.nii.gz -f
        scil_compute_mean_fixel_lobe_mrds_metric_from_bundles.py \$bundle \
            $pdd $md \${bname}_md_metric.nii.gz -f
        scil_compute_mean_fixel_lobe_mrds_metric_from_bundles.py \$bundle \
            $pdd $rd \${bname}_rd_metric.nii.gz -f
        scil_compute_mean_fixel_lobe_mrds_metric_from_bundles.py \$bundle \
            $pdd $ad \${bname}_ad_metric.nii.gz -f
    done
    """
}


process Bundle_Centroid {
    input:
    set sid, file(bundles) from bundles_for_centroids

    output:
    set sid, "*_centroid_${params.nb_points}.trk" into centroids_computed

    when:
    !params.use_provided_centroids

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    for bundle in $bundles_list;
        do if [[ \$bundle == *"__"* ]]; then
            pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
            bname=\${bundle:\$pos}
            bname=\$(basename \$bname .trk)
        else
            bname=\$(basename \$bundle .trk)
        fi
        bname=\${bname/_ic/}
        scil_compute_centroid.py \$bundle centroid.trk --nb_points $params.nb_points -f
        scil_uniformize_streamlines_endpoints.py centroid.trk ${sid}__\${bname}_centroid_${params.nb_points}.trk --auto
    done
    """
}

process Resample_Centroid {
    input:
    set sid, file(bundles) from centroids_for_resample

    output:
    set sid, "${sid}__*_centroid_${params.nb_points}.trk" into centroids_provided

    when:
    params.use_provided_centroids

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    for bundle in $bundles_list;
        do if [[ \$bundle == *"__"* ]]; then
            pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
            bname=\${bundle:\$pos}
            bname=\$(basename \$bname .trk)
        else
            bname=\$(basename \$bundle .trk)
        fi
        bname=\${bname/_centroid/}
        bname=\${bname/_ic/}

        echo \$(scil_resample_streamlines.py \$bundle \
            "${sid}__\${bname}_centroid_${params.nb_points}.trk" \
            --nb_pts_per_streamline $params.nb_points -f)
    done
    """
}

if (params.use_provided_centroids) {
    centroids_provided
        .set{centroids_for_label_and_distance_map}
}
else {
    centroids_computed
        .set{centroids_for_label_and_distance_map}
}

bundles_for_label_and_distance_map
    .join(centroids_for_label_and_distance_map, by: 0)
    .set{bundles_centroids_for_label_and_distance_map}

process Bundle_Label_And_Distance_Maps {
    input:
    set sid, file(bundles), file(centroids) from\
        bundles_centroids_for_label_and_distance_map

    output:
    set sid, "${sid}__*_labels.nii.gz", "${sid}__*_distances.nii.gz" into\
        label_distance_maps_for_mean_std_per_point
    set sid, "${sid}__*_labels.trk" into bundles_for_uniformize
    file "${sid}__*_distances.trk"
    set sid, "${sid}__*_labels.nii.gz" into voxel_label_maps_for_volume,
                                            voxel_label_map_for_lesion_load

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    for bundle in $bundles_list;
        do if [[ \$bundle == *"__"* ]]; then
            pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
            bname=\${bundle:\$pos}
            bname=\$(basename \$bname .trk)
        else
            bname=\$(basename \$bundle .trk)
        fi
        bname=\${bname/_ic/}

        centroid=${sid}__\${bname}_centroid_${params.nb_points}.trk
        if [[ -f \${centroid} ]]; then
            scil_compute_bundle_voxel_label_map.py \$bundle \${centroid} tmp_out -f

            mv tmp_out/labels_map.nii.gz ${sid}__\${bname}_labels.nii.gz
            mv tmp_out/distance_map.nii.gz ${sid}__\${bname}_distances.nii.gz

            mv tmp_out/labels.trk ${sid}__\${bname}_labels.trk
            mv tmp_out/distance.trk ${sid}__\${bname}_distances.trk
        fi

    done
    """
}

process Uniformize_Bundle {
    input:
    set sid, file(bundles) from bundles_for_uniformize

    output:
    set sid, "*_uniformized.trk" into bundles_for_coloring,
             bundles_for_lesion_load,
             bundles_for_mean_std, bundles_for_endpoints_map,
             bundles_for_endpoints_metrics, bundles_for_centroid,
             bundles_for_volume, bundles_for_streamline_count,
             bundles_for_mean_std_per_point, bundles_for_length_stats

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    for bundle in $bundles_list; do
        # Uniformize the bundle orientation as well as simplifying
        # filename convention
        if [[ \$bundle == *"__"* ]]; then
            scil_uniformize_streamlines_endpoints.py \$bundle \
                \${bundle/_labels.trk/_uniformized.trk} --auto -f
        else
            scil_uniformize_streamlines_endpoints.py \$bundle \
                ${sid}__\${bundle/_labels.trk/_uniformized.trk} --auto -f
        fi
    done
    """
}

lesion_for_lesion_load
    .join(bundles_for_lesion_load, by: 0)
    .join(voxel_label_map_for_lesion_load, by: 0)
    .set{lesion_bundles_voxel_label_maps_for_lesion_load}

process Lesion_Load {
    input:
    set sid, file(lesion), file(bundles), file(label_maps) from\
        lesion_bundles_voxel_label_maps_for_lesion_load

    output:
    file "${sid}__lesion_load.json" into lesion_load_to_aggregate
    file "${sid}__lesion_load_per_point.json" into lesion_load_per_point_to_aggregate
    set sid, "${sid}__lesion_load_per_point.json" into lesion_load_per_point_for_plot
    file "${sid}__lesion_streamlines_stats.json"
    file "${sid}__*_lesion_map.nii.gz"
    file "${sid}__lesion_stats.json"

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    mkdir streamlines_stats/
    mkdir lesion_load/
    mkdir lesion_load_per_point/
    for bundle in $bundles_list;
        do if [[ \$bundle == *"__"* ]]; then
            pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
            bname=\${bundle:\$pos}
            bname=\$(basename \$bname .trk)
        else
            bname=\$(basename \$bundle .trk)
        fi
        bname=\${bname/_uniformized/}
        rm ${sid}__lesion_stats.json -f
        mv ${sid}__\${bname}_labels.nii.gz \$bname.nii.gz
        mv \$bundle \$bname.trk

        scil_analyse_lesions_load.py $lesion lesion_load_per_point/\$bname.json \
            --bundle_labels_map \$bname.nii.gz \
            --out_lesion_atlas "${sid}__\${bname}_lesion_map.nii.gz" \
            --min_lesion_vol $params.min_lesion_vol

        scil_analyse_lesions_load.py $lesion lesion_load/\$bname.json \
            --bundle \$bname.trk --out_lesion_stats ${sid}__lesion_stats.json \
            --out_streamlines_stats streamlines_stats/\$bname.json \
            --min_lesion_vol $params.min_lesion_vol
    done

    scil_merge_json.py ${sid}__lesion_stats.json ${sid}__lesion_stats.json \
        --remove_parent_key --add_parent_key ${sid} -f

    cd streamlines_stats
    scil_merge_json.py *.json ../${sid}__lesion_streamlines_stats.json \
        --add_parent_key ${sid}

    cd ../lesion_load
    scil_merge_json.py *.json ../${sid}__lesion_load.json \
        --add_parent_key ${sid}

    cd ../lesion_load_per_point
    scil_merge_json.py *.json ../${sid}__lesion_load_per_point.json \
        --add_parent_key ${sid}
    """
}

process Color_Bundle {
    input:
    set sid, file(bundles) from bundles_for_coloring

    output:
    file "*_colored.trk"

    script:
    def json_str = JsonOutput.toJson(params.colors)
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    echo '$json_str' >> colors.json
    scil_assign_uniform_color_to_tractograms.py $bundles_list --dict_colors colors.json
    """
}

process Bundle_Length_Stats {
    input:
    set sid, file(bundles) from bundles_for_length_stats

    output:
    file "${sid}__length_stats.json" into bundles_length_stats_to_aggregate

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    for bundle in $bundles_list;
        do if [[ \$bundle == *"__"* ]]; then
            pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
            bname=\${bundle:\$pos}
            bname=\$(basename \$bname .trk)
        else
            bname=\$(basename \$bundle .trk)
        fi
        bname=\${bname/_uniformized/}
        scil_compute_streamlines_length_stats.py \$bundle > \$bname.json
        done

        scil_merge_json.py *.json ${sid}__length_stats.json --add_parent_key ${sid} \
            --keep_separate
    """
}

process Bundle_Endpoints_Map {
    input:
    set sid, file(bundles) from bundles_for_endpoints_map

    output:
    file "${sid}__endpoints_map_raw.json" into endpoints_maps_to_aggregate
    set sid, "*_endpoints_map_head.nii.gz", "*_endpoints_map_tail.nii.gz" \
        into endpoints_maps_for_roi_stats

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    for bundle in $bundles_list;
        do if [[ \$bundle == *"__"* ]]; then
            pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
            bname=\${bundle:\$pos}
            bname=\$(basename \$bname .trk)
        else
            bname=\$(basename \$bundle .trk)
        fi
        bname=\${bname/_uniformized/}
        mv \$bundle \$bname.trk

        scil_compute_endpoints_map.py \$bname.trk \
            ${sid}__\${bname}_endpoints_map_head.nii.gz \
            ${sid}__\${bname}_endpoints_map_tail.nii.gz >\
            ${sid}__\${bname}_endpoints_map_raw.json
    done
    scil_merge_json.py *_endpoints_map_raw.json ${sid}__endpoints_map_raw.json \
        --no_list --add_parent_key ${sid}
    """
}

metrics_for_endpoints_roi_stats
    .mix(fixel_afd_for_endpoints_roi_stats)
    .mix(fixel_fa_for_endpoints_roi_stats)
    .mix(fixel_md_for_endpoints_roi_stats)
    .mix(fixel_rd_for_endpoints_roi_stats)
    .mix(fixel_ad_for_endpoints_roi_stats)
    .groupTuple(by: 0)
    .map{it -> [it[0], it[1..-1].flatten()]}
    .set{metrics_afd_for_endpoints_roi_stats}

metrics_afd_for_endpoints_roi_stats
    .combine(endpoints_maps_for_roi_stats, by: 0)
    .set{metrics_endpoints_for_roi_stats}

process Bundle_Metrics_Stats_In_Endpoints {
    input:
    set sid, file(metrics), file(endpoints_map_head), file(endpoints_map_tail)\
         from metrics_endpoints_for_roi_stats

    output:
    file "${sid}__endpoints_metric_stats.json" into endpoints_metric_stats_to_aggregate

    script:
    normalize_weights =\
        params.endpoints_metric_stats_normalize_weights ?\
            '--normalize_weights' : '--bin'
    String map_list = endpoints_map_head.join(", ").replace(',', '')
    """
    shopt -s extglob
    for map in $map_list;
        do if [[ \$map == *"__"* ]]; then
            pos=\$((\$(echo \$map | grep -b -o __ | cut -d: -f1)+2))
            bname=\${map:\$pos}
            bname=\$(basename \$bname .nii.gz)
        else
            bname=\$(basename \$map .nii.gz)
        fi
        bname=\${bname/_endpoints_map_head/}
        mv \$map \${bname}_head.nii.gz
        mv \${map/_head/_tail} \${bname}_tail.nii.gz

        b_metrics="$metrics"
        b_metrics=\$(echo \$b_metrics | tr ' ' '\n' | grep -v "_afd_metric" | tr '\n' ' ')
        if [[ -f \${bname}_ic_afd_metric.nii.gz ]];
        then
            mv \${bname}_ic_afd_metric.nii.gz afd_metric.nii.gz
            b_metrics+=" afd_metric.nii.gz"
        fi

        scil_compute_metrics_stats_in_ROI.py \${bname}_head.nii.gz $normalize_weights\
            --metrics \${b_metrics} > \${bname}_head.json
        scil_compute_metrics_stats_in_ROI.py \${bname}_tail.nii.gz $normalize_weights\
            --metrics \${b_metrics} > \${bname}_tail.json
    done

    scil_merge_json.py *_tail.json *_head.json ${sid}__endpoints_metric_stats.json \
        --no_list --add_parent_key ${sid}
    """
}

metrics_for_endpoints_metrics
    .mix(fixel_afd_for_endpoints_metrics)
    .mix(fixel_fa_for_endpoints_metrics)
    .mix(fixel_md_for_endpoints_metrics)
    .mix(fixel_rd_for_endpoints_metrics)
    .mix(fixel_ad_for_endpoints_metrics)
    .groupTuple(by: 0)
   .map{it -> [it[0], it[1..-1].flatten()]}
    .set{metrics_afd_for_endpoints_metrics}

bundles_for_endpoints_metrics
    .flatMap{ sid, bundles -> bundles.collect{[sid, it]} }
    .combine(metrics_afd_for_endpoints_metrics, by: 0)
    .set{metrics_bundles_for_endpoints_metrics}

process Bundle_Endpoints_Metrics {
    input:
    set sid, file(bundle), file(metrics) from metrics_bundles_for_endpoints_metrics

    output:
    file "*/*_endpoints_metric.nii.gz"

    when:
    !params.skip_projection_endpoints_metrics

    script:
    """
    shopt -s extglob
    bundle=$bundle
    if [[ \$bundle == *"__"* ]]; then
        pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
        bname=\${bundle:\$pos}
        bname=\$(basename \$bname .trk)
    else
        bname=\$(basename \$bundle .trk)
    fi
    bname=\${bname/_uniformized/}
    bname=\${bname/_ic/}
    mkdir \${bname}

    b_metrics="$metrics"
    b_metrics=\$(echo \$b_metrics | tr ' ' '\n' | grep -v "_afd_metric" | tr '\n' ' ')
    if [[ -f \${bname}_ic_afd_metric.nii.gz ]];
    then
        mv \${bname}_ic_afd_metric.nii.gz afd_metric.nii.gz
        b_metrics+=" afd_metric.nii.gz"
    fi

    scil_project_streamlines_to_map.py \$bundle \${bname} --in_metrics \${b_metrics} --from_wm
    cd \${bname}
    for i in *.nii.gz;
        do mv "\$i" "${sid}__\$i";
    done

    for i in *;
        do
        if [[ \$i == "${sid}__${sid}__"* ]];
        then
            mv - \$i "\${i/${sid}__${sid}__/${sid}__/}";
        fi
    done
    """
}

metrics_for_mean_std
    .mix(fixel_afd_for_mean_std)
    .mix(fixel_fa_for_mean_std)
    .mix(fixel_md_for_mean_std)
    .mix(fixel_rd_for_mean_std)
    .mix(fixel_ad_for_mean_std)
    .groupTuple(by: 0)
    .map{it -> [it[0], it[1..-1].flatten()]}
    .set{metrics_afd_for_mean_std}

metrics_afd_for_mean_std
    .combine(bundles_for_mean_std, by: 0)
    .set{metrics_bundles_for_mean_std}
process Bundle_Mean_Std {
    input:
    set sid, file(metrics), file(bundles) from metrics_bundles_for_mean_std

    output:
    file "${sid}__mean_std.json" into mean_std_to_aggregate

    script:
    density_weighting =\
        params.mean_std_density_weighting ? '--density_weighting' : ''
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    shopt -s extglob
    for bundle in $bundles_list;
        do if [[ \$bundle == *"__"* ]]; then
            pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
            bname=\${bundle:\$pos}
            bname=\$(basename \$bname .trk)
        else
            bname=\$(basename \$bundle .trk)
        fi
        bname=\${bname/_uniformized/}
        bname=\${bname/_ic/}
        mv \$bundle \$bname.trk

        b_metrics="$metrics"
        b_metrics=\$(echo \$b_metrics | tr ' ' '\n' | grep -v "_afd_metric" | tr '\n' ' ')
        if [[ -f \${bname}_ic_afd_metric.nii.gz ]];
        then
            mv \${bname}_ic_afd_metric.nii.gz afd_metric.nii.gz
            b_metrics+=" afd_metric.nii.gz"
        fi

        scil_compute_bundle_mean_std.py $density_weighting \$bname.trk \${b_metrics} >\
            \${bname}.json
    done
    scil_merge_json.py *.json ${sid}__mean_std.json --no_list --add_parent_key ${sid}
    """
}

process Bundle_Volume {
    input:
    set sid, file(bundles) from bundles_for_volume

    output:
    file "${sid}__volume.json" into volumes_to_aggregate

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    for bundle in $bundles_list;
        do if [[ \$bundle == *"__"* ]]; then
            pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
            bname=\${bundle:\$pos}
            bname=\$(basename \$bname .trk)
        else
            bname=\$(basename \$bundle .trk)
        fi
        bname=\${bname/_uniformized/}
        mv \$bundle \$bname.trk
        scil_compute_bundle_volume.py \$bname.trk > \${bname}.json
    done
    scil_merge_json.py *.json ${sid}__volume.json --no_list --add_parent_key ${sid}
    """
}

process Bundle_Streamline_Count {
    input:
    set sid, file(bundles) from bundles_for_streamline_count

    output:
    file "${sid}__streamline_count.json" into streamline_counts_to_aggregate

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    for bundle in $bundles_list;
        do if [[ \$bundle == *"__"* ]]; then
            pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
            bname=\${bundle:\$pos}
            bname=\$(basename \$bname .trk)
        else
            bname=\$(basename \$bundle .trk)
        fi
        bname=\${bname/_uniformized/}
        mv \$bundle \$bname.trk
        scil_count_streamlines.py \$bname.trk > \${bname}.json
    done
    scil_merge_json.py *.json ${sid}__streamline_count.json --no_list \
        --add_parent_key ${sid}
    """
}


process Bundle_Volume_Per_Label {
    input:
    set sid, file(voxel_label_maps) from voxel_label_maps_for_volume

    output:
    file "${sid}__volume_per_label.json" into volumes_per_label_to_aggregate

    script:
    String maps_list = voxel_label_maps.join(", ").replace(',', '')
    """
    for map in $maps_list;
        do if [[ \$map == *"__"* ]]; then
            pos=\$((\$(echo \$map | grep -b -o __ | cut -d: -f1)+2))
            bname=\${map:\$pos}
            bname=\$(basename \$bname .nii.gz)
        else
            bname=\$(basename \$map .nii.gz)
        fi
        bname=\${bname/_voxel_label_map/}

        scil_compute_bundle_volume_per_label.py \$map \$bname --sort_keys >\
            \${bname}.json
        done
    scil_merge_json.py *.json ${sid}__volume_per_label.json --no_list \
        --add_parent_key ${sid}
    """
}

metrics_for_mean_std_per_point
    .mix(fixel_afd_for_mean_std_per_point)
    .mix(fixel_fa_for_mean_std_per_point)
    .mix(fixel_md_for_mean_std_per_point)
    .mix(fixel_rd_for_mean_std_per_point)
    .mix(fixel_ad_for_mean_std_per_point)
    .groupTuple(by: 0)
    .map{it -> [it[0], it[1..-1].flatten()]}
    .set{metrics_afd_for_std_per_point}

metrics_afd_for_std_per_point
    .join(bundles_for_mean_std_per_point, by: 0)
    .join(label_distance_maps_for_mean_std_per_point, by: 0)
    .set{metrics_bundles_label_distance_maps_for_mean_std_per_point}

process Bundle_Mean_Std_Per_Point {
    input:
    set sid, file(metrics), file(bundles), file(label_maps), file(distance_maps) \
         from metrics_bundles_label_distance_maps_for_mean_std_per_point

    output:
    set sid, "${sid}__mean_std_per_point.json" into \
        mean_std_per_point_for_plot
    file "${sid}__mean_std_per_point.json" into \
        mean_std_per_point_to_aggregate

    script:
    density_weighting =\
        params.mean_std_per_point_density_weighting ? '--density_weighting' : ''
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    shopt -s extglob
    for bundle in $bundles_list;
        do if [[ \$bundle == *"__"* ]]; then
            pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
            bname=\${bundle:\$pos}
            bname=\$(basename \$bname .trk)
        else
            bname=\$(basename \$bundle .trk)
        fi
        bname=\${bname/_uniformized/}
        bname=\${bname/_ic/} 
        mv \$bundle \$bname.trk
        label_map=${sid}__\${bname}_labels.nii.gz
        distance_map=${sid}__\${bname}_distances.nii.gz

        b_metrics="$metrics"
        b_metrics=\$(echo \$b_metrics | tr ' ' '\n' | grep -v "_afd_metric" | tr '\n' ' ')
        if [[ -f \${bname}_ic_afd_metric.nii.gz ]];
        then
            mv \${bname}_ic_afd_metric.nii.gz afd_metric.nii.gz
            b_metrics+=" afd_metric.nii.gz"
        fi

        scil_compute_bundle_mean_std_per_point.py \$bname.trk \$label_map \
            \${b_metrics} --sort_keys $density_weighting > \$bname.json
        done
        scil_merge_json.py *.json ${sid}__mean_std_per_point.json --no_list \
            --add_parent_key ${sid}
    """
}

mean_std_per_point_for_plot

process Plot_Mean_Std_Per_Point {
    input:
    set sid, file(mean_std_per_point) from mean_std_per_point_for_plot

    output:
    set sid, "*.png"

    script:
    def json_str = JsonOutput.toJson(params.colors)
    """
    echo '$json_str' >> colors.json
    scil_plot_mean_std_per_point.py $mean_std_per_point tmp_dir/ --dict_colors \
        colors.json --nb_pts $params.nb_points
    mv tmp_dir/* ./
    """
}

process Plot_Lesions_Per_Point {
    input:
    set sid, file(lesion_per_point) from lesion_load_per_point_for_plot

    output:
    set sid, "*.png"

    script:
    def json_str = JsonOutput.toJson(params.colors)
    """
    echo '$json_str' >> colors.json
    scil_merge_json.py $lesion_per_point tmp.json --recursive --average_last_layer
    scil_plot_mean_std_per_point.py tmp.json tmp_dir/ --dict_colors \
        colors.json --nb_pts $params.nb_points
    mv tmp_dir/* ./
    """
}
lesion_load_to_aggregate
    .collect()
    .set{all_lesion_load_to_aggregate}

process Aggregate_All_Lesion_Load {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from all_lesion_load_to_aggregate

    output:
    file "lesion_load.json"
    file "lesion_load.xlsx"

    script:
    """
    scil_merge_json.py $jsons lesion_load.json --average_last_layer --recursive
    scil_harmonize_json.py lesion_load.json lesion_load.json -f -v --sort_keys
    scil_convert_json_to_xlsx.py lesion_load.json lesion_load.xlsx
    """
}

lesion_load_per_point_to_aggregate
    .collect()
    .set{all_lesion_load_per_point_to_aggregate}

process Aggregate_All_Lesion_Load_Per_Point {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from all_lesion_load_per_point_to_aggregate

    output:
    file "lesion_load_per_point.json" into population_lesion_load_per_point
    file "lesion_load_per_point.xlsx"

    script:
    String json_list = jsons.join(", ").replace(',', '')
    """
    for json in $json_list
        do scil_merge_json.py \$json \${json/.json/_avg.json} --remove_parent_key --recursive --average_last_layer
    done
    scil_merge_json.py *_avg.json lesion_load_per_point.json  \
        --recursive
    scil_harmonize_json.py lesion_load_per_point.json lesion_load_per_point.json -f -v --sort_keys
    scil_convert_json_to_xlsx.py lesion_load_per_point.json lesion_load_per_point.xlsx \
        --stats_over_population
    """
}

endpoints_maps_to_aggregate
    .collect()
    .set{all_aggregate_endspoints_map}

process Aggregate_All_Endpoints_Map {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from all_aggregate_endspoints_map

    output:
    file "endpoints_map.json"
    file "endpoints_map.xlsx"

    script:
    """
    scil_merge_json.py $jsons endpoints_map.json --no_list
    scil_harmonize_json.py endpoints_map.json endpoints_map.json -f -v --sort_keys
    scil_convert_json_to_xlsx.py endpoints_map.json endpoints_map.xlsx
    """
}

endpoints_metric_stats_to_aggregate
    .collect()
    .set{all_aggregate_all_endpoints_metric_stats}

process Aggregate_All_Endpoints_Metric_Stats {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from all_aggregate_all_endpoints_metric_stats

    output:
    file "endpoints_metric_stats.json"
    file "endpoints_metric_stats.xlsx"

    script:
    """
    scil_merge_json.py $jsons endpoints_metric_stats.json --no_list
    scil_harmonize_json.py endpoints_metric_stats.json endpoints_metric_stats.json -f -v --sort_keys
    scil_convert_json_to_xlsx.py endpoints_metric_stats.json endpoints_metric_stats.xlsx
    """
}

bundles_length_stats_to_aggregate
    .collect()
    .set{all_bundle_length_stats_to_aggretate}

process Aggregate_All_Bundle_Length_Stats {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from all_bundle_length_stats_to_aggretate

    output:
    file "length_stats.json"
    file "length_stats.xlsx"

    script:
    """
    scil_merge_json.py $jsons length_stats.json --no_list
    scil_harmonize_json.py length_stats.json length_stats.json -f -v --sort_keys
    scil_convert_json_to_xlsx.py length_stats.json length_stats.xlsx
    """
}

mean_std_to_aggregate
    .collect()
    .set{all_mean_std_to_aggregate}

process Aggregate_All_mean_std {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from all_mean_std_to_aggregate

    output:
    file "mean_std.json"
    file "mean_std.xlsx"

    script:
    """
    scil_merge_json.py $jsons mean_std.json --no_list
    scil_harmonize_json.py mean_std.json mean_std.json -f -v --sort_keys
    scil_convert_json_to_xlsx.py mean_std.json mean_std.xlsx
    """
}

volumes_to_aggregate
    .collect()
    .set{all_volumes_to_aggregate}

process Aggregate_All_Volume {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from all_volumes_to_aggregate

    output:
    file "volumes.json"
    file "volumes.xlsx"

    script:
    """
    scil_merge_json.py $jsons volumes.json --no_list
    scil_harmonize_json.py volumes.json volumes.json -f -v --sort_keys
    scil_convert_json_to_xlsx.py volumes.json volumes.xlsx
    """
}

streamline_counts_to_aggregate
    .collect()
    .set{all_streamline_counts_to_aggregate}

process Aggregate_All_Streamline_Count {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from all_streamline_counts_to_aggregate

    output:
    file "streamline_count.json"
    file "streamline_count.xlsx"

    script:
    """
    scil_merge_json.py $jsons streamline_count.json --no_list
    scil_harmonize_json.py streamline_count.json streamline_count.json -f -v --sort_keys
    scil_convert_json_to_xlsx.py streamline_count.json streamline_count.xlsx
    """
}

volumes_per_label_to_aggregate
    .collect()
    .set{all_volumes_per_label_to_aggregate}

process Aggregate_All_Volume_Per_Label {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from all_volumes_per_label_to_aggregate

    output:
    file "volume_per_label.json"
    file "volume_per_label.xlsx"

    script:
    """
    scil_merge_json.py $jsons volume_per_label.json --no_list
    scil_harmonize_json.py volume_per_label.json volume_per_label.json -f -v --sort_keys
    scil_convert_json_to_xlsx.py volume_per_label.json volume_per_label.xlsx
    """
}

mean_std_per_point_to_aggregate
    .collect()
    .set{all_mean_std_per_point_to_aggregate}

process Aggregate_All_Mean_Std_Per_Point {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from all_mean_std_per_point_to_aggregate

    output:
    file "mean_std_per_point.json" into population_mean_std_per_point
    file "mean_std_per_point.xlsx"

    script:
    String json_list = jsons.join(", ").replace(',', '')
    """
    for json in $json_list
        do scil_merge_json.py \$json \${json/.json/_avg.json} --remove_parent_key --recursive
    done
    scil_merge_json.py *_avg.json mean_std_per_point.json  \
        --recursive
    scil_harmonize_json.py mean_std_per_point.json mean_std_per_point.json -f -v --sort_keys
    scil_convert_json_to_xlsx.py mean_std_per_point.json mean_std_per_point.xlsx \
        --stats_over_population
    """
}

population_mean_std_per_point
    .concat(population_lesion_load_per_point)
    .set{population_mean_std_lesion_per_point}

process Plot_Population_Mean_Std_Per_Point {
    tag = { "Plots" }
    publishDir = params.plotPublishDir

    input:
    file(json_a) from population_mean_std_lesion_per_point

    output:
    file "*.png"

    script:
    def json_str = JsonOutput.toJson(params.colors)
    """
    echo '$json_str' >> colors.json
    scil_plot_mean_std_per_point.py $json_a tmp_dir/ --dict_colors colors.json \
        --stats_over_population --nb_pts $params.nb_points
    mv tmp_dir/* ./
    """
}
