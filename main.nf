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
                "mean_std_density_weighting":"$params.mean_std_density_weighting",
                "mean_std_per_point_density_weighting":"$params.mean_std_per_point_density_weighting",
                "endpoints_metric_stats_normalize_weights":"$params.endpoints_metric_stats_normalize_weights",
                "voxel_label_map_upsample":"$params.voxel_label_map_upsample",
                "cpu_count":"$cpu_count",
                "skip_pruning":"$params.skip_pruning",
                "skip_outlier_rejection":"$params.skip_outlier_rejection"
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
    .into{bundles_for_coloring;bundles_for_centroids}

Channel
    .fromFilePairs("$params.input/**/metrics/*.nii.gz",
                   size: -1) { it.parent.parent.name }
    .set{in_metrics}

Channel
    .fromFilePairs("$params.input/**/centroids/*.trk",
        size: -1) { it.parent.parent.name }
    .into{centroids_for_resample; in_centroids_check}

in_metrics.into{metrics_for_mean_std;
                metrics_for_endpoints_metrics; metrics_for_endpoints_roi_stats;
                metrics_for_volume;
                metrics_for_mean_std_per_point}

process Bundle_Centroid {
    input:
    set sid, file(bundles) from bundles_for_centroids

    output:
    set sid, "*_centroid.trk" into centroids_computed

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
    
        scil_compute_centroid.py \$bundle centroid.trk --nb_points $params.nb_points -f
        scil_uniformize_streamlines_endpoints.py centroid.trk ${sid}__\${bname}_centroid.trk --auto
    done
    """
}

process Resample_Centroid {
    input:
    set sid, file(bundles) from centroids_for_resample

    output:
    set sid, "${sid}__*_centroid.trk" into centroids_provided

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
    
        scil_resample_streamlines.py \$bundle "${sid}__\${bname}_centroid.trk" --nb_pts_per_streamline $params.nb_points
    done
    """
}

if (params.use_provided_centroids) {
    centroids_provided
        .into{centroids_for_label_and_distance_map; centroids_for_volume_per_label;lol}
}
else {
    centroids_computed
        .into{centroids_for_label_and_distance_map; centroids_for_volume_per_label;lol}
}

process Color_Bundle {
    input:
    set sid, file(bundles) from bundles_for_coloring

    output:
    set sid, "*_colored.trk" into bundles_for_mean_std, bundles_for_endpoints_map,
             bundles_for_endpoints_metrics, bundles_for_centroid,
             bundles_for_label_and_distance_map, bundles_for_volume,
             bundles_for_streamline_count, bundles_for_volume_per_label,
             bundles_for_mean_std_per_point, bundles_for_length_stats

    script:
    def json_str = JsonOutput.toJson(params.colors)
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    echo '$json_str' >> colors.json
    scil_assign_color_to_trk.py $bundles_list --dict_colors colors.json
    for bundle in *_colored.trk; do
        if [[ \$bundle == *"__"* ]]; then
            scil_uniformize_streamlines_endpoints.py \$bundle \$bundle --auto -f
        else
            scil_uniformize_streamlines_endpoints.py \$bundle ${sid}__\$bundle --auto -f
            rm \$bundle
        fi
    done
    """
}

process Bundle_Length_Stats {
    input:
    set sid, file(bundles) from bundles_for_length_stats

    output:
    set sid, "${sid}__length_stats.json" into bundle_length_stats_to_aggregate

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
        bname=\${bname/_colored/}
        scil_compute_streamlines_length_stats.py \$bundle > \$bname.json
        done

        scil_merge_json.py *.json ${sid}__length_stats.json --add_parent_key ${sid} --keep_separate
    """
}

process Bundle_Endpoints_Map {
    input:
    set sid, file(bundles) from bundles_for_endpoints_map

    output:
    set sid, "${sid}__endpoints_map_raw.json" into endpoints_map_to_aggregate
    set sid, "*_endpoints_map_head.nii.gz", "*_endpoints_map_tail.nii.gz" \
        into endpoints_map_for_roi_stats

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
        bname=\${bname/_colored/}
        mv \$bundle \$bname.trk

        scil_compute_endpoints_map.py \$bname.trk ${sid}__\${bname}_endpoints_map_head.nii.gz \
            ${sid}__\${bname}_endpoints_map_tail.nii.gz >\
                ${sid}__\${bname}_endpoints_map_raw.json
    done
    scil_merge_json.py *_endpoints_map_raw.json ${sid}__endpoints_map_raw.json --no_list --add_parent_key ${sid}
    """
}
metrics_for_endpoints_roi_stats
    .combine(endpoints_map_for_roi_stats, by: 0)
    .set{metrics_endpoints_for_roi_stats}

process Bundle_Metrics_Stats_In_Endpoints {
    input:
    set sid, file(metrics), file(endpoints_map_head), file(endpoints_map_tail) \
         from metrics_endpoints_for_roi_stats

    output:
    set sid, "${sid}__endpoints_metric_stats.json" into\
        endpoints_metric_stats_to_aggregate

    script:
    normalize_weights =\
        params.endpoints_metric_stats_normalize_weights ?\
            '--normalize_weights' : '--bin'
    String map_list = endpoints_map_head.join(", ").replace(',', '')
    """
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

        scil_compute_metrics_stats_in_ROI.py \${bname}_head.nii.gz $normalize_weights\
            --metrics $metrics > \${bname}_head.json
        scil_compute_metrics_stats_in_ROI.py \${bname}_tail.nii.gz $normalize_weights\
            --metrics $metrics > \${bname}_tail.json
    done

    scil_merge_json.py *_tail.json *_head.json ${sid}__endpoints_metric_stats.json --no_list --add_parent_key ${sid}
    """
}

metrics_for_endpoints_metrics
    .combine(bundles_for_endpoints_metrics, by: 0)
    .set{metrics_bundles_for_endpoints_metrics}

process Bundle_Endpoints_Metrics {
    input:
    set sid, file(metrics), file(bundles) from metrics_bundles_for_endpoints_metrics

    output:
    file "*/*_endpoints_metric.nii.gz"

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
        bname=\${bname/_colored/}
        mkdir \${bname}

        scil_compute_endpoints_metric.py \$bundle $metrics \${bname}
        cd \${bname}
        for i in *.nii.gz; do mv "\$i" "${sid}__\$i"; done
        rename s/${sid}__${sid}__/${sid}__/ *
        cd ../
    done
    """
}

metrics_for_mean_std
    .combine(bundles_for_mean_std, by: 0)
    .set{metrics_bundles_for_mean_std}
process Bundle_Mean_Std {
    input:
    set sid, file(metrics), file(bundles) from metrics_bundles_for_mean_std

    output:
    set sid, "${sid}__mean_std.json" into mean_std_to_aggregate

    script:
    density_weighting =\
        params.mean_std_density_weighting ? '--density_weighting' : ''
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
        bname=\${bname/_colored/}
        mv \$bundle \$bname.trk
        scil_compute_bundle_mean_std.py $density_weighting \$bname.trk $metrics >\
            \${bname}.json
    done
    scil_merge_json.py *.json ${sid}__mean_std.json --no_list --add_parent_key ${sid}
    """
}

bundles_for_label_and_distance_map
    .join(centroids_for_label_and_distance_map, by: 0)
    .set{bundles_centroids_for_label_and_distance_map}

process Bundle_Label_And_Distance_Maps {
    input:
    set sid, file(bundles), file(centroids) from\
        bundles_centroids_for_label_and_distance_map

    output:
    set sid, "${sid}__*_labels.npz", "${sid}__*_distances.npz" into\
        label_distance_maps_for_mean_std_per_point

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
        bname=\${bname/_colored/}
        centroid=\${bundle/_colored/_centroid}
        scil_label_and_distance_maps.py \$bundle \$centroid\
            ${sid}__\${bname}_labels.npz ${sid}__\${bname}_distances.npz
        done
    """
}

process Bundle_Volume {
    input:
    set sid, file(bundles) from bundles_for_volume

    output:
    set sid, "${sid}__volume.json" into volume_to_aggregate

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
        bname=\${bname/_colored/}
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
    set sid, "${sid}__streamline_count.json" into streamline_count_to_aggregate

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
        bname=\${bname/_colored/}
        mv \$bundle \$bname.trk
        scil_count_streamlines.py \$bname.trk > \${bname}.json
    done
    scil_merge_json.py *.json ${sid}__streamline_count.json --no_list --add_parent_key ${sid}
    """
}

bundles_for_volume_per_label
    .join(centroids_for_volume_per_label)
    .set{bundles_centroids_voxel_label_map}

process Bundle_Voxel_Label_Map {
    input:
    set sid, file(bundles), file(centroid) from\
        bundles_centroids_voxel_label_map

    output:
    set sid, "${sid}__*_voxel_label_map.nii.gz" into\
        voxel_label_maps_for_volume

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
        bname=\${bname/_colored/}
        centroid=\${bundle/_colored/_centroid}
        scil_compute_bundle_voxel_label_map.py  \$bundle \$centroid \
            ${sid}__\${bname}_voxel_label_map.nii.gz
    done
    """
}

process Bundle_Volume_Per_Label {
    input:
    set sid, file(voxel_label_maps) from voxel_label_maps_for_volume

    output:
    set sid, "${sid}__volume_per_label.json" into volume_per_label_to_aggregate

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
    scil_merge_json.py *.json ${sid}__volume_per_label.json --no_list --add_parent_key ${sid}
    """
}

metrics_for_mean_std_per_point
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
    for bundle in $bundles_list;
        do if [[ \$bundle == *"__"* ]]; then
            pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
            bname=\${bundle:\$pos}
            bname=\$(basename \$bname .trk)
        else
            bname=\$(basename \$bundle .trk)
        fi
        bname=\${bname/_colored/}
        mv \$bundle \$bname.trk
        label_map=${sid}__\${bname}_labels.npz
        distance_map=${sid}__\${bname}_labels.npz

        scil_compute_bundle_mean_std_per_point.py \$bname.trk \$label_map \$distance_map \
            $metrics --sort_keys $density_weighting > \$bname.json
        done
        scil_merge_json.py *.json ${sid}__mean_std_per_point.json --no_list --add_parent_key ${sid}
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
    scil_plot_mean_std_per_point.py $mean_std_per_point tmp_dir/ --dict_colors colors.json
    mv tmp_dir/* ./
    """
}

// bundle_length_stats_to_aggregate
//     .groupTuple()
//     .set{sorted_bundle_length_stats_to_aggretate}

// process Aggregate_Subject_Bundle_Length_Stats {
//     tag = { "${sid}" }

//     input:
//     set sid, file(jsons) from sorted_bundle_length_stats_to_aggretate

//     output:
//     file "${sid}__bundle_length_stats.json"\
//         into all_bundle_length_stats_to_aggretate

//     script:
//     """
//     cat *.json | jq -s 'reduce .[] as \$item ({}; . * \$item)' >\
//         ${sid}__bundle_length_stats.json
//     """
// }

// endpoints_map_to_aggregate
//     .groupTuple()
//     .map{ch1, ch2, ch3 -> [ch1, ch3]}
//     .set{sorted_endpoints_map_to_aggregate}

// process Aggregate_Subject_Endpoints_Map {
//     tag = { "${sid}" }

//     input:
//     set sid, file(jsons) from sorted_endpoints_map_to_aggregate

//     output:
//     file "${sid}__endpoints_map_agg.json" into all_endpoints_map_to_aggregate

//     script:
//     """
//     cat *endpoints_map.json | jq -s '[.[] | to_entries] | flatten |\
//         reduce .[] as \$dot ({}; .[\$dot.key] += \$dot.value)' >\
//             ${sid}__endpoints_map_agg.json
//     """
// }

// endpoints_metric_stats_to_aggregate
//     .groupTuple()
//     .map{ch1, ch2, ch3 -> [ch1, ch3]}
//     .set{sorted_endpoints_metric_stats_to_aggregate}

// process Aggregate_Subject_Endpoints_Metric_Stats {
//     tag = { "${sid}" }

//     input:
//     set sid, file(jsons) from sorted_endpoints_metric_stats_to_aggregate

//     output:
//     file "${sid}__endpoints_metric_stats_agg.json" into\
//         all_endpoints_metric_stats_to_aggregate

//     script:
//     """
//     cat *endpoints_metric_stats.json | jq -s '[.[] | to_entries] | flatten |\
//         reduce .[] as \$dot ({}; .[\$dot.key] += \$dot.value)' >\
//             ${sid}__endpoints_metric_stats_agg.json
//     """
// }

// mean_std_to_aggregate
//     .groupTuple()
//     .map{ch1, ch2, ch3 -> [ch1, ch3]}
//     .set{sorted_mean_std_to_aggregate}

// process Aggregate_Subject_mean_std {
//     tag = { "${sid}" }

//     input:
//     set sid, file(jsons) from sorted_mean_std_to_aggregate

//     output:
//     file "${sid}__mean_std_agg.json" into all_mean_std_to_aggregate

//     script:
//     """
//     cat *mean_std.json | jq -s '[.[] | to_entries] | flatten |\
//         reduce .[] as \$dot ({}; .[\$dot.key] += \$dot.value)' >\
//             ${sid}__mean_std_agg.json
//     """
// }

// volume_to_aggregate
//     .groupTuple()
//     .map{ch1, ch2, ch3 -> [ch1, ch3]}
//     .set{sorted_volume_to_aggregate}

// process Aggregate_Subject_Volume {
//     tag = { "${sid}" }

//     input:
//     set sid, file(jsons) from sorted_volume_to_aggregate

//     output:
//     file "${sid}__volume_agg.json" into all_volume_to_aggregate

//     script:
//     """
//     cat *volume.json | jq -s '[.[] | to_entries] | flatten |\
//         reduce .[] as \$dot ({}; .[\$dot.key] += \$dot.value)' >\
//             ${sid}__volume_agg.json
//     """
// }

// streamline_count_to_aggregate
//     .groupTuple()
//     .map{ch1, ch2, ch3 -> [ch1, ch3]}
//     .set{sorted_streamline_count_to_aggregate}

// process Aggregate_Subject_streamline_count {
//     tag = { "${sid}" }

//     input:
//     set sid, file(jsons) from sorted_streamline_count_to_aggregate

//     output:
//     file "${sid}__streamline_count_agg.json" into all_streamline_count_to_aggregate

//     script:
//     """
//     cat *streamline_count.json | jq -s '[.[] | to_entries] | flatten |\
//         reduce .[] as \$dot ({}; .[\$dot.key] += \$dot.value)' >\
//             ${sid}__streamline_count_agg.json
//     """
// }

// volume_per_label_to_aggregate
//     .groupTuple()
//     .map{ch1, ch2, ch3 -> [ch1, ch3]}
//     .set{sorted_volume_per_label_to_aggregate}

// process Aggregate_Subject_Volume_Per_Label {
//     tag = { "${sid}" }

//     input:
//     set sid, file(jsons) from sorted_volume_per_label_to_aggregate

//     output:
//     file "${sid}__volume_per_label_agg.json" into\
//         all_volume_per_label_to_aggregate

//     script:
//     """
//     cat *volume_per_label.json | jq -s '[.[] | to_entries] | flatten |\
//         reduce .[] as \$dot ({}; .[\$dot.key] += \$dot.value)' >\
//             ${sid}__volume_per_label_agg.json
//     """
// }

// all_endpoints_map_to_aggregate
//     .collect()
//     .set{in_aggregate_all_endpoints_map}

// process Aggregate_All_Endpoints_Map {
//     tag = { "Statistics" }
//     publishDir = params.statsPublishDir

//     input:
//     file jsons from in_aggregate_all_endpoints_map

//     output:
//     file "endpoints_map.json"

//     """
//     cat *.json | jq -s add > endpoints_map.json
//     """
// }

// all_endpoints_metric_stats_to_aggregate
//     .collect()
//     .set{in_aggregate_all_endpoints_metric_stats}

// process Aggregate_All_Endpoints_Metric_Stats {
//     tag = { "Statistics" }
//     publishDir = params.statsPublishDir

//     input:
//     file jsons from in_aggregate_all_endpoints_metric_stats

//     output:
//     file "endpoints_metric_stats.json"

//     script:
//     """
//     cat *.json | jq -s add > endpoints_metric_stats.json
//     """
// }

// all_bundle_length_stats_to_aggretate
//     .collect()
//     .set{in_aggregate_all_bundle_length_stats}

// process Aggregate_All_Bundle_Length_Stats {
//     tag = { "Statistics" }
//     publishDir = params.statsPublishDir

//     input:
//     file jsons from in_aggregate_all_bundle_length_stats

//     output:
//     file "length_stats.json"

//     script:
//     """
//     cat *.json | jq -s add > length_stats.json
//     """
// }

// all_mean_std_to_aggregate
//     .collect()
//     .set{in_aggregate_all_mean_std}

// process Aggregate_All_mean_std {
//     tag = { "Statistics" }
//     publishDir = params.statsPublishDir

//     input:
//     file jsons from in_aggregate_all_mean_std

//     output:
//     file "mean_std.json"

//     """
//     cat *.json | jq -s add > mean_std.json
//     """
// }

// all_volume_to_aggregate
//     .collect()
//     .set{in_aggregate_all_volume}

// process Aggregate_All_Volume {
//     tag = { "Statistics" }
//     publishDir = params.statsPublishDir

//     input:
//     file jsons from in_aggregate_all_volume

//     output:
//     file "volumes.json"

//     """
//     cat *.json | jq -s add > volumes.json
//     """
// }

// all_streamline_count_to_aggregate
//     .collect()
//     .set{in_aggregate_all_streamline_count}

// process Aggregate_All_streamline_count {
//     tag = { "Statistics" }
//     publishDir = params.statsPublishDir

//     input:
//     file jsons from in_aggregate_all_streamline_count

//     output:
//     file "streamline_count.json"

//     """
//     cat *.json | jq -s add > streamline_count.json
//     """
// }

// all_volume_per_label_to_aggregate
//     .collect()
//     .set{in_aggregate_all_volume_per_label}

// process Aggregate_All_Volume_Per_Label {
//     tag = { "Statistics" }
//     publishDir = params.statsPublishDir

//     input:
//     file jsons from in_aggregate_all_volume_per_label

//     output:
//     file "volume_per_label.json"

//     """
//     cat *.json | jq -s add > volume_per_label.json
//     """
// }

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

    """
    scil_merge_json.py $jsons "mean_std_per_point.json" --remove_parent_key --recursive
    """
}

process Plot_Population_Mean_Std_Per_Point {
    tag = { "Plots" }
    publishDir = params.plotPublishDir

    input:
    file(json) from population_mean_std_per_point

    output:
    file "*.png"

    script:
    def json_str = JsonOutput.toJson(params.colors)
    """
    echo '$json_str' >> colors.json
    scil_plot_mean_std_per_point.py $json tmp_dir/ --dict_colors colors.json --stats_over_population
    mv tmp_dir/* ./
    """
}