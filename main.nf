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
    .into{bundles_for_uniformize;bundles_for_centroids}

Channel
    .fromFilePairs("$params.input/**/metrics/*.nii.gz",
                   size: -1) { it.parent.parent.name }
    .set{in_metrics}

Channel
    .fromFilePairs("$params.input/**/centroids/*.trk",
        size: -1) { it.parent.parent.name }
    .into{centroids_for_resample; in_centroids_check}

in_metrics
    .set{metrics_for_rename}

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

process Bundle_Centroid {
    input:
    set sid, file(bundles) from bundles_for_centroids

    output:
    set sid, "*_centroid_${$params.nb_points}.trk" into centroids_computed

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
        bname=\${bname/$params.bundle_suffix_to_remove/}
    
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

        scil_resample_streamlines.py \$bundle \
            "${sid}__\${bname}_centroid_${params.nb_points}.trk" \
            --nb_pts_per_streamline $params.nb_points -f
    done
    """
}

if (params.use_provided_centroids) {
    centroids_provided
        .into{centroids_for_label_and_distance_map;centroids_for_volume_per_label}
}
else {
    centroids_computed
        .into{centroids_for_label_and_distance_map;centroids_for_volume_per_label}
}

process Uniformize_Bundle {
    input:
    set sid, file(bundles) from bundles_for_uniformize

    output:
    set sid, "*_uniformized.trk" into bundles_for_coloring
    set sid, "*_uniformized.trk" into bundles_for_mean_std, bundles_for_endpoints_map,
             bundles_for_endpoints_metrics, bundles_for_centroid,
             bundles_for_label_and_distance_map, bundles_for_volume,
             bundles_for_streamline_count, bundles_for_volume_per_label,
             bundles_for_mean_std_per_point, bundles_for_length_stats

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    for bundle in $bundles_list; do
        # Uniformize the bundle orientation as well as simplifying
        # filename convention
        if [[ \$bundle == *"__"* ]]; then
            scil_uniformize_streamlines_endpoints.py \$bundle \
                \${bundle/.trk/_uniformized.trk} --auto -f
        else
            scil_uniformize_streamlines_endpoints.py \$bundle \
                ${sid}__\${bundle/.trk/_uniformized.trk} --auto -f
        fi
    done

    # Remove suffix from RecobundlesX if present
    for bundle in *_uniformized.trk; do
        mv \$bundle \${bundle/$params.bundle_suffix_to_remove/}
    done
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
    scil_assign_color_to_trk.py $bundles_list --dict_colors colors.json
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
    .combine(endpoints_maps_for_roi_stats, by: 0)
    .set{metrics_endpoints_for_roi_stats}

process Bundle_Metrics_Stats_In_Endpoints {
    input:
    set sid, file(metrics), file(endpoints_map_head), file(endpoints_map_tail) \
         from metrics_endpoints_for_roi_stats

    output:
    file "${sid}__endpoints_metric_stats.json" into endpoints_metric_stats_to_aggregate

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

    scil_merge_json.py *_tail.json *_head.json ${sid}__endpoints_metric_stats.json \
        --no_list --add_parent_key ${sid}
    """
}

bundles_for_endpoints_metrics
    .flatMap{ sid, bundles -> bundles.collect{[sid, it]} }
    .combine(metrics_for_endpoints_metrics, by: 0)
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
    bundle=$bundle
    if [[ \$bundle == *"__"* ]]; then
        pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
        bname=\${bundle:\$pos}
        bname=\$(basename \$bname .trk)
    else
        bname=\$(basename \$bundle .trk)
    fi
    bname=\${bname/_uniformized/}
    mkdir \${bname}

    scil_compute_endpoints_metric.py \$bundle $metrics \${bname}
    cd \${bname}
    for i in *.nii.gz; do mv "\$i" "${sid}__\$i"; done
    rename s/${sid}__${sid}__/${sid}__/ *

    """
}

metrics_for_mean_std
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
        bname=\${bname/_uniformized/}
        centroid=\${bundle/_uniformized/_centroid_${params.nb_points}}
        scil_label_and_distance_maps.py \$bundle \$centroid\
            ${sid}__\${bname}_labels.npz ${sid}__\${bname}_distances.npz
        done
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

bundles_for_volume_per_label
    .join(centroids_for_volume_per_label)
    .set{bundles_centroids_voxel_label_map}

process Bundle_Voxel_Label_Map {
    input:
    set sid, file(bundles), file(centroid) from\
        bundles_centroids_voxel_label_map

    output:
    set sid, "${sid}__*_voxel_label_map.nii.gz" into voxel_label_maps_for_volume

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
        centroid=\${bundle/_uniformized/_centroid_${params.nb_points}}
        scil_compute_bundle_voxel_label_map.py  \$bundle \$centroid \
            ${sid}__\${bname}_voxel_label_map.nii.gz
    done
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
        bname=\${bname/_uniformized/}
        mv \$bundle \$bname.trk
        label_map=${sid}__\${bname}_labels.npz
        distance_map=${sid}__\${bname}_labels.npz

        scil_compute_bundle_mean_std_per_point.py \$bname.trk \$label_map \$distance_map \
            $metrics --sort_keys $density_weighting > \$bname.json
        done
        scil_merge_json.py *.json ${sid}__mean_std_per_point.json --no_list \
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
    scil_plot_mean_std_per_point.py $mean_std_per_point tmp_dir/ --dict_colors \
        colors.json
    mv tmp_dir/* ./
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

    """
    scil_merge_json.py $jsons endpoints_map.json --no_list
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

    """
    scil_merge_json.py $jsons mean_std.json --no_list
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

    """
    scil_merge_json.py $jsons volumes.json --no_list
    scil_convert_json_to_xlsx.py volumes.json volumes.xlsx
    """
}

streamline_counts_to_aggregate
    .collect()
    .set{all_streamline_counts_to_aggregate}

process Aggregate_All_streamline_count {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from all_streamline_counts_to_aggregate

    output:
    file "streamline_count.json"
    file "streamline_count.xlsx"

    """
    scil_merge_json.py $jsons streamline_count.json --no_list
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

    """
    scil_merge_json.py $jsons volume_per_label.json --no_list
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

    """
    scil_merge_json.py $jsons "mean_std_per_point.json" --remove_parent_key \
        --recursive
    scil_convert_json_to_xlsx.py mean_std_per_point.json mean_std_per_point.xlsx \
        --stats_over_population
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
    scil_plot_mean_std_per_point.py $json tmp_dir/ --dict_colors colors.json \
        --stats_over_population
    mv tmp_dir/* ./
    """
}