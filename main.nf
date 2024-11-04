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
                "skip_tract_profiles": "$params.skip_tract_profiles",
                "cpu_count":"$cpu_count",
                "bundle_suffix_to_remove":"$params.bundle_suffix_to_remove"
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

in_metrics
    .set{metrics_for_rename}


in_bundles_check.map{it[1]}.flatten().count().set{number_bundles_for_compare}
in_centroids_check.map{it[1]}.flatten().count().set{number_centroids_for_compare}

if (params.use_provided_centroids && !params.skip_tract_profiles){
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
    set sid, "${sid}__*_ic.trk" into bundles_for_label_and_distance_map, bundles_for_centroids, bundles_for_fixel_afd, bundles_if_skip_tract_profiles

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

      scil_tractogram_remove_invalid.py \$bundle ${sid}__\${bname}_ic.trk --remove_single_point --remove_overlapping_points --cut_invalid --no_empty
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
        scil_bundle_mean_fixel_afd.py \$bundle $fodf \${bname}_afd_metric.nii.gz
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
    !params.skip_tract_profiles

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
        scil_bundle_compute_centroid.py \$bundle centroid.trk --nb_points $params.nb_points -f
        scil_bundle_uniformize_endpoints.py centroid.trk ${sid}__\${bname}_centroid_${params.nb_points}.trk --auto
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
    !params.skip_tract_profiles

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

        echo \$(scil_tractogram_resample_nb_points.py \$bundle \
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
    set sid, "${sid}__*_labels.trk" into bundles_labels_for_uniformize
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
            scil_bundle_label_map.py \$bundle \${centroid} tmp_out -f

            mv tmp_out/labels_map.nii.gz ${sid}__\${bname}_labels.nii.gz
            mv tmp_out/distance_map.nii.gz ${sid}__\${bname}_distances.nii.gz

            mv tmp_out/labels.trk ${sid}__\${bname}_labels.trk
            mv tmp_out/distance.trk ${sid}__\${bname}_distances.trk
        fi
    done
    """
}

if (params.skip_tract_profiles) {
    bundles_if_skip_tract_profiles
        .set{bundles_for_uniformize}
}
else {
    bundles_labels_for_uniformize
        .set{bundles_for_uniformize}
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
            if [[ \$bundle == *"labels"* ]]; then
                scil_bundle_uniformize_endpoints.py \$bundle \
                    \${bundle/_labels.trk/_uniformized.trk} --auto -f
            elif [[ \$bundle == *"ic"* ]]; then
                scil_bundle_uniformize_endpoints.py \$bundle \
                    \${bundle/_ic.trk/_uniformized.trk} --auto -f
            fi
        else
            if [[ \$bundle == *"labels"* ]]; then
                scil_bundle_uniformize_endpoints.py \$bundle \
                    ${sid}__\${bundle/_labels.trk/_uniformized.trk} --auto -f
            elif [[ \$bundle == *"ic"* ]]; then
                scil_bundle_uniformize_endpoints.py \$bundle \
                    ${sid}__\${bundle/_ic.trk/_uniformized.trk} --auto -f
            fi
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
    scil_labels_from_mask.py $lesion lesion_labels.nii.gz
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

        scil_lesions_info.py lesion_labels.nii.gz lesion_load_per_point/\$bname.json \
            --bundle_labels_map \$bname.nii.gz \
            --out_lesion_atlas "${sid}__\${bname}_lesion_map.nii.gz" \
            --min_lesion_vol $params.min_lesion_vol

        scil_lesions_info.py lesion_labels.nii.gz lesion_load/\$bname.json \
            --bundle \$bname.trk --out_lesion_stats ${sid}__lesion_stats.json \
            --out_streamlines_stats streamlines_stats/\$bname.json \
            --min_lesion_vol $params.min_lesion_vol
    done

    scil_json_merge_entries.py ${sid}__lesion_stats.json ${sid}__lesion_stats.json \
        --remove_parent_key --add_parent_key ${sid} -f

    cd streamlines_stats
    scil_json_merge_entries.py *.json ../${sid}__lesion_streamlines_stats.json \
        --add_parent_key ${sid}

    cd ../lesion_load
    scil_json_merge_entries.py *.json ../${sid}__lesion_load.json \
        --add_parent_key ${sid}

    cd ../lesion_load_per_point
    scil_json_merge_entries.py *.json ../${sid}__lesion_load_per_point.json \
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
    scil_tractogram_assign_uniform_color.py $bundles_list --dict_colors colors.json --out_suffix colored
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
        scil_tractogram_print_info.py \$bundle > \$bname.json
        done

        scil_json_merge_entries.py *.json ${sid}__length_stats.json --add_parent_key ${sid} \
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

        scil_bundle_compute_endpoints_map.py \$bname.trk \
            ${sid}__\${bname}_endpoints_map_head.nii.gz \
            ${sid}__\${bname}_endpoints_map_tail.nii.gz >\
            ${sid}__\${bname}_endpoints_map_raw.json
    done
    scil_json_merge_entries.py *_endpoints_map_raw.json ${sid}__endpoints_map_raw.json \
        --no_list --add_parent_key ${sid}
    """
}

metrics_for_endpoints_roi_stats
    .mix(fixel_afd_for_endpoints_roi_stats)
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

        scil_volume_stats_in_ROI.py \${bname}_head.nii.gz $normalize_weights\
            --metrics \${b_metrics} > \${bname}_head.json
        scil_volume_stats_in_ROI.py \${bname}_tail.nii.gz $normalize_weights\
            --metrics \${b_metrics} > \${bname}_tail.json
    done

    scil_json_merge_entries.py *_tail.json *_head.json ${sid}__endpoints_metric_stats.json \
        --no_list --add_parent_key ${sid}
    """
}

metrics_for_endpoints_metrics
    .mix(fixel_afd_for_endpoints_metrics)
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

    scil_tractogram_project_streamlines_to_map.py \$bundle \${bname} --in_metrics \${b_metrics} --from_wm
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

        scil_bundle_mean_std.py $density_weighting \$bname.trk \${b_metrics} >\
            \${bname}.json
    done
    scil_json_merge_entries.py *.json ${sid}__mean_std.json --no_list --add_parent_key ${sid}
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
        scil_bundle_shape_measures.py \$bname.trk --out_json \${bname}.json
    done
    scil_json_merge_entries.py *.json ${sid}__volume.json --no_list --add_parent_key ${sid} --keep_separate
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

        scil_bundle_volume_per_label.py \$map \$bname --sort_keys >\
            \${bname}.json
        done
    scil_json_merge_entries.py *.json ${sid}__volume_per_label.json --no_list \
        --add_parent_key ${sid}
    """
}

metrics_for_mean_std_per_point
    .mix(fixel_afd_for_mean_std_per_point)
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

        scil_bundle_mean_std.py \$bname.trk --per_point \$label_map \
            \${b_metrics} --sort_keys $density_weighting > \$bname.json
        done
        scil_json_merge_entries.py *.json ${sid}__mean_std_per_point.json --no_list \
            --add_parent_key ${sid}
    """
}

process Plot_Mean_Std_Per_Point {
    input:
    set sid, file(mean_std_per_point) from mean_std_per_point_for_plot

    output:
    set sid, "*.png"

    script:
    def json_str = JsonOutput.toJson(params.colors)
    """
    echo '$json_str' >> colors.json
    scil_plot_stats_per_point.py $mean_std_per_point tmp_dir/ --dict_colors \
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
    scil_json_merge_entries.py $lesion_per_point tmp.json --recursive --average_last_layer
    scil_plot_stats_per_point.py tmp.json tmp_dir/ --dict_colors \
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
    scil_json_merge_entries.py $jsons lesion_load.json --average_last_layer --recursive
    scil_json_harmonize_entries.py lesion_load.json lesion_load.json -f -v --sort_keys
    scil_json_convert_entries_to_xlsx.py lesion_load.json lesion_load.xlsx
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
        do scil_json_merge_entries.py \$json \${json/.json/_avg.json} --remove_parent_key --recursive --average_last_layer
    done
    scil_json_merge_entries.py *_avg.json lesion_load_per_point.json  \
        --recursive
    scil_json_harmonize_entries.py lesion_load_per_point.json lesion_load_per_point.json -f -v --sort_keys
    scil_json_convert_entries_to_xlsx.py lesion_load_per_point.json lesion_load_per_point.xlsx \
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
    scil_json_merge_entries.py $jsons endpoints_map.json --no_list
    scil_json_harmonize_entries.py endpoints_map.json endpoints_map.json -f -v --sort_keys
    scil_json_convert_entries_to_xlsx.py endpoints_map.json endpoints_map.xlsx
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
    scil_json_merge_entries.py $jsons endpoints_metric_stats.json --no_list
    scil_json_harmonize_entries.py endpoints_metric_stats.json endpoints_metric_stats.json -f -v --sort_keys
    scil_json_convert_entries_to_xlsx.py endpoints_metric_stats.json endpoints_metric_stats.xlsx
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
    scil_json_merge_entries.py $jsons length_stats.json --no_list
    scil_json_harmonize_entries.py length_stats.json length_stats.json -f -v --sort_keys
    scil_json_convert_entries_to_xlsx.py length_stats.json length_stats.xlsx
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
    scil_json_merge_entries.py $jsons mean_std.json --no_list
    scil_json_harmonize_entries.py mean_std.json mean_std.json -f -v --sort_keys
    scil_json_convert_entries_to_xlsx.py mean_std.json mean_std.xlsx
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
    scil_json_merge_entries.py $jsons volumes.json --no_list
    scil_json_harmonize_entries.py volumes.json volumes.json -f -v --sort_keys
    scil_json_convert_entries_to_xlsx.py volumes.json volumes.xlsx
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
    scil_json_merge_entries.py $jsons volume_per_label.json --no_list
    scil_json_harmonize_entries.py volume_per_label.json volume_per_label.json -f -v --sort_keys
    scil_json_convert_entries_to_xlsx.py volume_per_label.json volume_per_label.xlsx
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
        do scil_json_merge_entries.py \$json \${json/.json/_avg.json} --remove_parent_key --recursive
    done
    scil_json_merge_entries.py *_avg.json mean_std_per_point.json  \
        --recursive
    scil_json_harmonize_entries.py mean_std_per_point.json mean_std_per_point.json -f -v --sort_keys
    scil_json_convert_entries_to_xlsx.py mean_std_per_point.json mean_std_per_point.xlsx \
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
    scil_plot_stats_per_point.py $json_a tmp_dir/ --dict_colors colors.json \
        --stats_over_population --nb_pts $params.nb_points
    mv tmp_dir/* ./
    """
}
