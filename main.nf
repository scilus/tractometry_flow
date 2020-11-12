#!/usr/bin/env nextflow
import com.google.common.hash.Hashing
import com.google.common.base.Charsets

params.help = false
params.root = false

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
    .fromFilePairs("$params.root/**/bundles/*.trk",
                   size: -1) { it.parent.parent.name }
    .set{in_bundles}

Channel
    .fromFilePairs("$params.root/**/metrics/*.nii.gz",
                   size: -1) { it.parent.parent.name }
    .set{in_metrics}

Channel
    .fromFilePairs("$params.root/**/centroids/*.trk",
        size: -1) { it.parent.parent.name }
    .into{in_centroids; in_centroids_check}

in_metrics.into{metrics_for_mean_std;
                metrics_for_endpoints_metrics; metrics_for_endpoints_roi_stats;
                metrics_for_volume;
                metrics_for_mean_std_per_point}

in_bundles
    .flatMap{ sid, bundles -> bundles.collect{ [sid, it] } }
    .set{bundles_to_prune}

in_centroids
    .flatMap{ sid, centroid -> centroid.collect {
            [sid, it.name.find(~/.*(?=_centroid\.trk)/), it] } }
    .set{centroids_for_resample}

process Prune_Bundle {
    input:
    set sid, file(bundle) from bundles_to_prune

    output:
    set sid, val(bname), "${sid}__${bname}__pruned.trk" optional true into\
        pruned_bundle_for_outliers

    script:
    bname = bundle.name.take(bundle.name.lastIndexOf('.'))
    if (bname.contains('__'))
    {
        bname = bname.substring(bname.lastIndexOf("__") + 2)
    }
    min_length = params.pruning.default.min_length
    max_length = params.pruning.default.max_length
    for (entry in params.pruning) {
        def match = bundle.find { it =~ entry.key }
        if (match) {
            min_length = entry.value.min_length
            max_length = entry.value.max_length
            break
        }
    }

    if (!params.skip_pruning) {
        """
        scil_filter_streamlines_by_length.py --minL $min_length --maxL $max_length \
            $bundle ${sid}__${bname}__pruned.trk
        """
    }
    else {
        """
        cp -Lr $bundle ${sid}__${bname}__pruned.trk
        """
    }
}

process Remove_Bundle_Outliers {
    input:
    set sid, val(bname), file(bundle) from pruned_bundle_for_outliers

    output:
    set sid, val(bname), "${sid}__${bname}__outliers_removed.trk" optional true\
        into outlier_rejected_bundle_for_color
    file "${sid}__${bname}__outliers.trk" optional true

    script:
    alpha = params.outlier_rejection.default.alpha
    for (entry in params.outlier_rejection) {
        def match = bundle.find { it =~ entry.key }
        if (match) {
            alpha = entry.value.alpha
            break
        }
    }

    if (!params.skip_outlier_rejection) {
        """
        scil_outlier_rejection.py $bundle ${sid}__${bname}__outliers_removed.trk \
            --remaining_bundle ${sid}__${bname}__outliers.trk --alpha $alpha
        """
    } else {
        """
        cp -Lr $bundle ${sid}__${bname}__outliers_removed.trk
        """
    }
}

process Color_Bundle {
    input:
    set sid, val(bname), file(bundle) from outlier_rejected_bundle_for_color

    output:
    set sid, val(bname), "${sid}__${bname}__outliers_removed_colored.trk" optional true\
        into outlier_rejected_bundle_for_mean_std,
             outlier_rejected_bundle_for_endpoints_map,
             outlier_rejected_bundle_for_endpoints_metrics,
             outlier_rejected_bundle_for_centroid,
             outlier_rejected_bundle_for_label_and_distance_map,
             outlier_rejected_bundle_for_volume,
             outlier_rejected_bundle_for_streamline_count,
             outlier_rejected_bundle_for_volume_per_label,
             outlier_rejected_bundle_for_mean_std_per_point,
             outlier_rejected_bundle_for_length_stats
    set sid, val(bname), val(color) optional true into color_for_mean_std_plot
    set val(bname), val(color) optional true into color_for_population_mean_std_plot

    script:
    Random rand_color_gen = new Random(
        Hashing.md5().hashString(bname, Charsets.UTF_8).asLong())
    double rnd_num = rand_color_gen.nextDouble()
    color = hsv_to_rgb(rnd_num % 1, 0.99, 0.99)
    for (entry in params.colors) {
        def match = bundle.find { it =~ entry.key }
        if (match) {
            color = entry.value
            break
        }
    }
    """
    scil_assign_color_to_trk.py $bundle\
        ${sid}__${bname}__outliers_removed_colored.trk\
        $color
    """
}

process Bundle_Length_Stats {
    input:
    set sid, bname, file(bundle) from outlier_rejected_bundle_for_length_stats

    output:
    set sid, "${sid}__${bname}__length_stats.json" into\
        bundle_length_stats_to_aggregate

    script:
    """
    scil_compute_streamlines_length_stats.py $bundle >\
        ${sid}__${bname}__length_stats_raw.json
    jq '{"${sid}": {"${bname}": .}}' \
        ${sid}__${bname}__length_stats_raw.json >\
            ${sid}__${bname}__length_stats.json
    """
}

process Bundle_Endpoints_Map {
    input:
    set sid, bname, file(bundle) from outlier_rejected_bundle_for_endpoints_map

    output:
    set sid, val(bname), "${sid}__${bname}__endpoints_map.json"\
        into endpoints_map_to_aggregate
    set sid, val(bname), "${sid}__${bname}__endpoints_map_head.nii.gz",\
        "${sid}__${bname}__endpoints_map_tail.nii.gz"\
        into endpoints_map_for_roi_stats

    script:
    """
    scil_endpoints_map.py $bundle ${sid}__${bname}__endpoints_map_head.nii.gz \
        ${sid}__${bname}__endpoints_map_tail.nii.gz >\
            ${sid}__${bname}__endpoints_map_raw.json
    jq 'to_entries|map({(.key|ltrimstr("${sid}__")|\
                         sub("__outliers_removed_colored"; "")):.value})|\
        reduce .[] as \$item ({}; . * \$item)|\
        {"${sid}": .}' ${sid}__${bname}__endpoints_map_raw.json >\
            ${sid}__${bname}__endpoints_map.json
    """
}

metrics_for_endpoints_roi_stats
    .cross(endpoints_map_for_roi_stats)
    .map{ch1, ch2 -> [*ch2, *ch1[1..-1]]}
    .set{in_endpoints_roi_stats}

process Bundle_Metrics_Stats_In_Endpoints {
    input:
    set sid, val(bname), file(endpoints_map_head), file(endpoints_map_tail),
        file(metrics) from in_endpoints_roi_stats

    output:
    set sid, val(bname), "${sid}__${bname}__endpoints_metric_stats.json" into\
        endpoints_metric_stats_to_aggregate

    script:
    normalize_weights =\
        params.endpoints_metric_stats_normalize_weights ?\
            '--normalize_weights' : '--bin'
    """
    scil_compute_metrics_stats_in_ROI.py $endpoints_map_head $normalize_weights\
        --metrics $metrics\
        > ${sid}__${bname}__endpoints_head_metric_stats_raw.json
    scil_compute_metrics_stats_in_ROI.py $endpoints_map_tail $normalize_weights\
        --metrics $metrics\
        > ${sid}__${bname}__endpoints_tail_metric_stats_raw.json
    jq 'to_entries|map({(.key|ltrimstr("${sid}__")|\
                         sub("__endpoints_map"; "")):.value})|\
        reduce .[] as \$item ({}; . + \$item)|\
        {"${sid}": .}' ${sid}__${bname}__endpoints_head_metric_stats_raw.json\
            ${sid}__${bname}__endpoints_tail_metric_stats_raw.json >\
            ${sid}__${bname}__endpoints_metric_stats_raw.json
    cat ${sid}__${bname}__endpoints_metric_stats_raw.json |\
        jq -s '[.[] | to_entries] | flatten |\
            reduce .[] as \$dot ({}; .[\$dot.key] += \$dot.value)' >\
        ${sid}__${bname}__endpoints_metric_stats.json
    """
}

metrics_for_endpoints_metrics
    .cross(outlier_rejected_bundle_for_endpoints_metrics)
    .map{ch1, ch2 -> [*ch2, *ch1[1..-1]]}
    .set{in_endpoints_metrics}

process Bundle_Endpoints_Metrics {
    input:
    set sid, bname, file(bundle), file(metrics) from in_endpoints_metrics

    output:
    file "${bname}/${sid}__*endpoints_metric.nii.gz"

    script:
    """
    mkdir ${bname}
    scil_endpoints_metric.py $bundle $metrics ${bname}
    cd ${bname}
    for i in *.nii.gz; do mv "\$i" "${sid}__\$i"; done
    rename s/${sid}__${sid}__/${sid}__/ *
    """
}

metrics_for_mean_std
    .cross(outlier_rejected_bundle_for_mean_std)
    .map{ch1, ch2 -> [*ch2, *ch1[1..-1]]}
    .set{in_bundle_mean_std}

process Bundle_Mean_Std {
    input:
    set sid, bname, file(bundle), file(metrics) from in_bundle_mean_std

    output:
    set sid, val(bname), "${sid}__${bname}__mean_std.json" into\
        mean_std_to_aggregate

    script:
    density_weighting =\
        params.mean_std_density_weighting ? '--density_weighting' : ''
    """
    scil_bundle_mean_std.py $density_weighting $bundle $metrics >\
        ${sid}__${bname}__mean_std_raw.json
    jq 'to_entries|map({(.key|ltrimstr("${sid}__")|\
                         rtrimstr("__outliers_removed_colored")):.value})|\
        reduce .[] as \$item ({}; . + \$item)|\
        {"${sid}": .}' ${sid}__${bname}__mean_std_raw.json >\
            ${sid}__${bname}__mean_std.json

    """
}

process Bundle_Centroid {
    input:
    set sid, bname, file(bundle) from outlier_rejected_bundle_for_centroid

    output:
    set sid, val(bname), "${sid}__${bname}__centroid.trk" into centroids_computed

    when:
    !params.use_provided_centroids

    script:
    """
    scil_compute_centroid.py $bundle centroid.trk --nb_points $params.nb_points
    scil_uniformize_streamlines_endpoints.py centroid.trk ${sid}__${bname}__centroid.trk --auto
    """
}

process Resample_Centroid {
    input:
    set sid, bname, file(bundle) from centroids_for_resample

    output:
    set sid, val(bname), "${sid}__${bname}__centroid.trk" into centroids_provided

    when:
    params.use_provided_centroids

    script:
    bname = bname.replace("_centroid", "")
    if (bname.contains('__'))
    {
        bname = bname.substring(bname.lastIndexOf("__") + 2)
    }
    """
    scil_resample_streamlines.py $bundle "${sid}__${bname}__centroid.trk" --nb_pts_per_streamline $params.nb_points
    """
}

if (params.use_provided_centroids) {
    centroids_provided
        .into{centroid_for_label_and_distance_map; centroid_for_volume_per_label;lol}
}
else {
    centroids_computed
        .into{centroid_for_label_and_distance_map; centroid_for_volume_per_label;lol}
}

outlier_rejected_bundle_for_label_and_distance_map
    .phase(centroid_for_label_and_distance_map)
        {it -> it[0] + it[1] }
    .map{ch1, ch2 -> [*ch1, ch2[2]]}
    .set{bundle_outlier_rejected_and_centroid_for_label_and_distance_map}

process Bundle_Label_And_Distance_Maps {
    input:
    set sid, bname, file(bundle), file(centroid) from\
        bundle_outlier_rejected_and_centroid_for_label_and_distance_map

    output:
    set sid, val(bname), "${sid}__${bname}__labels.npz",
        "${sid}__${bname}__distances.npz" into\
        label_distance_maps_for_mean_std_per_point

    script:
    """
    scil_label_and_distance_maps.py $bundle $centroid\
        ${sid}__${bname}__labels.npz ${sid}__${bname}__distances.npz
    """
}

process Bundle_Volume {
    input:
    set sid, bname, file(bundle) from outlier_rejected_bundle_for_volume

    output:
    set sid, val(bname), "${sid}__${bname}__volume.json" into\
        volume_to_aggregate

    script:
    """
    scil_compute_bundle_volume.py $bundle > ${sid}__${bname}__volume_raw.json
    jq 'to_entries|map({(.key|ltrimstr("${sid}__")|\
                         rtrimstr("__outliers_removed_colored")):.value})|\
        reduce .[] as \$item ({}; . + \$item)|\
        {"${sid}": .}' ${sid}__${bname}__volume_raw.json >\
            ${sid}__${bname}__volume.json
    """
}

process Bundle_streamline_count {
    input:
    set sid, bname, file(bundle) from outlier_rejected_bundle_for_streamline_count

    output:
    set sid, val(bname), "${sid}__${bname}__streamline_count.json" into\
        streamline_count_to_aggregate

    script:
    """
    scil_count_streamlines.py $bundle > ${sid}__${bname}__streamline_count_raw.json
    jq 'to_entries|map({(.key|ltrimstr("${sid}__")|\
                         rtrimstr("__outliers_removed_colored")):.value})|\
        reduce .[] as \$item ({}; . + \$item)|\
        {"${sid}": .}' ${sid}__${bname}__streamline_count_raw.json >\
            ${sid}__${bname}__streamline_count.json
    """
}

outlier_rejected_bundle_for_volume_per_label
    .phase(centroid_for_volume_per_label)
        {it -> it[0] + it[1]}
    .map{ch1, ch2 -> [*ch1, ch2[2]]}
    .set{in_bundle_voxel_label_map}

process Bundle_Voxel_Label_Map {
    input:
    set sid, bname, file(bundle), file(centroid) from\
        in_bundle_voxel_label_map

    output:
    set sid, val(bname), "${sid}__${bname}__voxel_label_map.nii.gz" into\
        voxel_label_map_for_volume

    script:
    """
    scil_bundle_voxel_label_map.py  $bundle $centroid ${sid}__${bname}__voxel_label_map.nii.gz
    """
}

process Bundle_Volume_Per_Label {
    input:
    set sid, bname, file(voxel_label_map) from voxel_label_map_for_volume

    output:
    set sid, val(bname), "${sid}__${bname}__volume_per_label.json" into\
            volume_per_label_to_aggregate

    script:
    """
    scil_bundle_volume_per_label.py\
        --sort_keys\
        $voxel_label_map $bname >\
        ${sid}__${bname}__volume_per_label_raw.json
    jq 'to_entries|map({(.key|ltrimstr("${sid}__")|\
                         rtrimstr("__outliers_removed_colored")):.value})|\
        reduce .[] as \$item ({}; . + \$item)|\
        {"${sid}": .}' ${sid}__${bname}__volume_per_label_raw.json >\
            ${sid}__${bname}__volume_per_label.json
    """
}

outlier_rejected_bundle_for_mean_std_per_point
    .phase(label_distance_maps_for_mean_std_per_point)
        {it -> it[0] + it[1]}
    .map{ch1, ch2 -> [*ch1, *ch2[2..3]]}
    .set{bundle_and_label_distance_maps_for_mean_std_per_point}

metrics_for_mean_std_per_point
    .cross(bundle_and_label_distance_maps_for_mean_std_per_point)
    .map{ch1, ch2 -> [*ch2, *ch1[1..-1]]}
    .set{in_bundle_mean_std_per_point}

process Bundle_Mean_Std_Per_Point {
    input:
    set sid, bname, file(bundle), file(label_map), file(distance_map),
        file(metrics) from in_bundle_mean_std_per_point

    output:
    set sid, val(bname), "${sid}__${bname}__mean_std_per_point_raw.json" into\
        mean_std_per_point_for_plot
    set sid, val(bname), "${sid}__${bname}__mean_std_per_point.json" into\
        mean_std_per_point_to_aggregate

    script:
    density_weighting =\
        params.mean_std_per_point_density_weighting ? '--density_weighting' : ''
    """
    scil_bundle_mean_std_per_point.py --sort_keys $density_weighting $bundle\
        $label_map $distance_map $metrics |\
        jq 'to_entries|map({(.key|rtrimstr("__outliers_removed_colored")):.value})|
                           reduce .[] as \$item ({}; . + \$item)'\
        > ${sid}__${bname}__mean_std_per_point_raw.json
    jq 'to_entries|map({(.key|ltrimstr("${sid}__")):.value})|\
        reduce .[] as \$item ({}; . + \$item)|\
        {"${sid}": .}' ${sid}__${bname}__mean_std_per_point_raw.json >\
            ${sid}__${bname}__mean_std_per_point.json
    """
}

mean_std_per_point_for_plot
    .combine(color_for_mean_std_plot, by: [0, 1])
    .set{mean_std_per_point_and_color_for_plot}

process Bundle_Plot_Mean_Std_Per_Point {
    input:
    set sid, val(bname), file(mean_std_per_point),
        val(color) from mean_std_per_point_and_color_for_plot

    output:
    set sid, val(bname), "*.png"

    script:
    """
    scil_plot_mean_std_per_point.py $mean_std_per_point tmp_dir/ --fill_color $color
    mv tmp_dir/* ./
    """
}

bundle_length_stats_to_aggregate
    .groupTuple()
    .set{sorted_bundle_length_stats_to_aggretate}

process Aggregate_Subject_Bundle_Length_Stats {
    tag = { "${sid}" }

    input:
    set sid, file(jsons) from sorted_bundle_length_stats_to_aggretate

    output:
    file "${sid}__bundle_length_stats.json"\
        into all_bundle_length_stats_to_aggretate

    script:
    """
    cat *.json | jq -s 'reduce .[] as \$item ({}; . * \$item)' >\
        ${sid}__bundle_length_stats.json
    """
}

endpoints_map_to_aggregate
    .groupTuple()
    .map{ch1, ch2, ch3 -> [ch1, ch3]}
    .set{sorted_endpoints_map_to_aggregate}

process Aggregate_Subject_Endpoints_Map {
    tag = { "${sid}" }

    input:
    set sid, file(jsons) from sorted_endpoints_map_to_aggregate

    output:
    file "${sid}__endpoints_map_agg.json" into all_endpoints_map_to_aggregate

    script:
    """
    cat *endpoints_map.json | jq -s '[.[] | to_entries] | flatten |\
        reduce .[] as \$dot ({}; .[\$dot.key] += \$dot.value)' >\
            ${sid}__endpoints_map_agg.json
    """
}

endpoints_metric_stats_to_aggregate
    .groupTuple()
    .map{ch1, ch2, ch3 -> [ch1, ch3]}
    .set{sorted_endpoints_metric_stats_to_aggregate}

process Aggregate_Subject_Endpoints_Metric_Stats {
    tag = { "${sid}" }

    input:
    set sid, file(jsons) from sorted_endpoints_metric_stats_to_aggregate

    output:
    file "${sid}__endpoints_metric_stats_agg.json" into\
        all_endpoints_metric_stats_to_aggregate

    script:
    """
    cat *endpoints_metric_stats.json | jq -s '[.[] | to_entries] | flatten |\
        reduce .[] as \$dot ({}; .[\$dot.key] += \$dot.value)' >\
            ${sid}__endpoints_metric_stats_agg.json
    """
}

mean_std_to_aggregate
    .groupTuple()
    .map{ch1, ch2, ch3 -> [ch1, ch3]}
    .set{sorted_mean_std_to_aggregate}

process Aggregate_Subject_mean_std {
    tag = { "${sid}" }

    input:
    set sid, file(jsons) from sorted_mean_std_to_aggregate

    output:
    file "${sid}__mean_std_agg.json" into all_mean_std_to_aggregate

    script:
    """
    cat *mean_std.json | jq -s '[.[] | to_entries] | flatten |\
        reduce .[] as \$dot ({}; .[\$dot.key] += \$dot.value)' >\
            ${sid}__mean_std_agg.json
    """
}

volume_to_aggregate
    .groupTuple()
    .map{ch1, ch2, ch3 -> [ch1, ch3]}
    .set{sorted_volume_to_aggregate}

process Aggregate_Subject_Volume {
    tag = { "${sid}" }

    input:
    set sid, file(jsons) from sorted_volume_to_aggregate

    output:
    file "${sid}__volume_agg.json" into all_volume_to_aggregate

    script:
    """
    cat *volume.json | jq -s '[.[] | to_entries] | flatten |\
        reduce .[] as \$dot ({}; .[\$dot.key] += \$dot.value)' >\
            ${sid}__volume_agg.json
    """
}

streamline_count_to_aggregate
    .groupTuple()
    .map{ch1, ch2, ch3 -> [ch1, ch3]}
    .set{sorted_streamline_count_to_aggregate}

process Aggregate_Subject_streamline_count {
    tag = { "${sid}" }

    input:
    set sid, file(jsons) from sorted_streamline_count_to_aggregate

    output:
    file "${sid}__streamline_count_agg.json" into all_streamline_count_to_aggregate

    script:
    """
    cat *streamline_count.json | jq -s '[.[] | to_entries] | flatten |\
        reduce .[] as \$dot ({}; .[\$dot.key] += \$dot.value)' >\
            ${sid}__streamline_count_agg.json
    """
}

volume_per_label_to_aggregate
    .groupTuple()
    .map{ch1, ch2, ch3 -> [ch1, ch3]}
    .set{sorted_volume_per_label_to_aggregate}

process Aggregate_Subject_Volume_Per_Label {
    tag = { "${sid}" }

    input:
    set sid, file(jsons) from sorted_volume_per_label_to_aggregate

    output:
    file "${sid}__volume_per_label_agg.json" into\
        all_volume_per_label_to_aggregate

    script:
    """
    cat *volume_per_label.json | jq -s '[.[] | to_entries] | flatten |\
        reduce .[] as \$dot ({}; .[\$dot.key] += \$dot.value)' >\
            ${sid}__volume_per_label_agg.json
    """
}

mean_std_per_point_to_aggregate
    .groupTuple()
    .map{ch1, ch2, ch3 -> [ch1, ch3]}
    .set{sorted_mean_std_per_point_to_aggregate}

process Aggregate_Subject_Mean_Std_Per_Point {
    tag = { "${sid}" }

    input:
    set sid, file(jsons) from sorted_mean_std_per_point_to_aggregate

    output:
    file "${sid}__mean_std_per_point_agg.json" into\
        all_mean_std_per_point_to_aggregate

    script:
    """
    cat *mean_std_per_point.json | jq -s '[.[] | to_entries] | flatten |\
        reduce .[] as \$dot ({}; .[\$dot.key] += \$dot.value)' >\
            ${sid}__mean_std_per_point_agg.json
    """
}

all_endpoints_map_to_aggregate
    .collect()
    .set{in_aggregate_all_endpoints_map}

process Aggregate_All_Endpoints_Map {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from in_aggregate_all_endpoints_map

    output:
    file "endpoints_map.json"

    """
    cat *.json | jq -s add > endpoints_map.json
    """
}

all_endpoints_metric_stats_to_aggregate
    .collect()
    .set{in_aggregate_all_endpoints_metric_stats}

process Aggregate_All_Endpoints_Metric_Stats {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from in_aggregate_all_endpoints_metric_stats

    output:
    file "endpoints_metric_stats.json"

    script:
    """
    cat *.json | jq -s add > endpoints_metric_stats.json
    """
}

all_bundle_length_stats_to_aggretate
    .collect()
    .set{in_aggregate_all_bundle_length_stats}

process Aggregate_All_Bundle_Length_Stats {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from in_aggregate_all_bundle_length_stats

    output:
    file "length_stats.json"

    script:
    """
    cat *.json | jq -s add > length_stats.json
    """
}

all_mean_std_to_aggregate
    .collect()
    .set{in_aggregate_all_mean_std}

process Aggregate_All_mean_std {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from in_aggregate_all_mean_std

    output:
    file "mean_std.json"

    """
    cat *.json | jq -s add > mean_std.json
    """
}

all_volume_to_aggregate
    .collect()
    .set{in_aggregate_all_volume}

process Aggregate_All_Volume {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from in_aggregate_all_volume

    output:
    file "volumes.json"

    """
    cat *.json | jq -s add > volumes.json
    """
}

all_streamline_count_to_aggregate
    .collect()
    .set{in_aggregate_all_streamline_count}

process Aggregate_All_streamline_count {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from in_aggregate_all_streamline_count

    output:
    file "streamline_count.json"

    """
    cat *.json | jq -s add > streamline_count.json
    """
}

all_volume_per_label_to_aggregate
    .collect()
    .set{in_aggregate_all_volume_per_label}

process Aggregate_All_Volume_Per_Label {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from in_aggregate_all_volume_per_label

    output:
    file "volume_per_label.json"

    """
    cat *.json | jq -s add > volume_per_label.json
    """
}

all_mean_std_per_point_to_aggregate
    .collect()
    .set{in_aggregate_all_mean_std_per_point}

process Aggregate_All_Mean_Std_Per_Point {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file jsons from in_aggregate_all_mean_std_per_point

    output:
    file "mean_std_per_point.json" into in_population_mean_std_per_point

    """
    cat *.json | jq -s add > mean_std_per_point.json
    """
}

process Population_Mean_Std_Per_Point {
    tag = { "Statistics" }
    publishDir = params.statsPublishDir

    input:
    file json from in_population_mean_std_per_point

    output:
    file "population_mean_std_per_point.json" into population_mean_std_per_point

    script:
    """
    jq -S 'def merge(b): (objects | reduce (b|to_entries[]) as \$e\
            (.; .[\$e.key] |= merge(\$e.value))) // .+b;\
        .[][][][][]|=[.] | reduce .[] as \$v ({}; merge(\$v)) |\
        .[][][][] |= add/length' $json > population_mean_std_per_point.json
    """
}

color_for_population_mean_std_plot
    .unique()
    .collate(2)
    .collect()
    .set{unique_color_for_population_mean_std_plot}

population_mean_std_per_point
    .combine(unique_color_for_population_mean_std_plot)
    .map{ch1 -> [ch1[0], ch1[1..-1]]}
    .set{population_mean_std_per_point_and_colors}

process Plot_Population_Mean_Std_Per_Point {
    tag = { "Plots" }
    publishDir = params.plotPublishDir + "/Mean_Std_Per_Point"

    input:
    set file(json), val(colors) from population_mean_std_per_point_and_colors

    output:
    file "*.png"

    script:
    def colormap = colors.collectEntries()
    bash_map = "( "
    for (e in colormap) {
        bash_map += "[" + "\"" + e.key + "\"]=\"" + e.value + "\" "
    }
    bash_map += ")"
    """
    declare -A lut=${bash_map}
    mkdir out_dir/
    for i in \$(jq -r 'keys[]' $json);
    do
        jq -r --arg k \$i '[{key: \$k, value: .[\$k]}] |\
            from_entries' $json > tmp_\${i}.json
        scil_plot_mean_std_per_point.py tmp_\${i}.json tmp_dir/ --fill_color \${lut[\$i]} -f
        mv tmp_dir/* ./
    done
    """
}

// ***************************
// **** UTILITY FUNCTIONS ****
// ***************************

def rgbToString(r, g, b) {
    String rs = String.format("%2s",
        Integer.toHexString((int)(r * 256))).replace(' ', '0');
    String gs = String.format("%2s",
        Integer.toHexString((int)(g * 256))).replace(' ', '0');
    String bs = String.format("%2s",
        Integer.toHexString((int)(b * 256))).replace(' ', '0');
    return "0x" + rs + gs + bs;
}

def hsv_to_rgb(hue, saturation, value) {
    h = (int)(hue * 6)
    float f = hue * 6 - h
    float p = value * (1 - saturation)
    float q = value * (1 - f * saturation)
    float t = value * (1 - (1 - f) * saturation)

    switch (h) {
      case 0: return rgbToString(value, t, p)
      case 1: return rgbToString(q, value, p)
      case 2: return rgbToString(p, value, t)
      case 3: return rgbToString(p, q, value)
      case 4: return rgbToString(t, p, value)
      case 5: return rgbToString(value, p, q)
      default: throw new RuntimeException(
        "Something went wrong when converting from HSV to RGB")
    }
}
