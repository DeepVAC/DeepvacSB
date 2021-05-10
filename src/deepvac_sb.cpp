/*
 * Copyright (c) 2020 DeepVAC
 * This file is part of libdeepvac, licensed under the GPLv3 (the "License")
 * You may not use this file except in compliance with the License.
 */
#include <cmath>
#include "deepvac_sb.h"

constexpr int DeepvacSB::resolution_ratio_level_[]{3200, 2100, 1700, 1200, 900, 600, 400};
constexpr int DeepvacSB::downsample_factor_[]{12, 8, 6, 5, 4, 3, 2, 1};

std::vector<unsigned int> DeepvacSB::getSbdIdx() noexcept {
    std::vector<unsigned int> sbd_indexs;
    //auto start = std::chrono::system_clock::now();
    auto size = hist_info_.size();
    if(!size){
        return sbd_indexs;
    }
    auto count = size / DeepvacSB::hist_size_ / DeepvacSB::channels_;
    if(count == 1) {
        sbd_indexs.push_back(0);
        return sbd_indexs;
    }
    std::vector<unsigned int>  hist_key_indexs = filterHistFrameIndex();
    std::vector<unsigned int>  hsv_key_indexs = filterHsvFrameIndex();

    sbd_indexs = getSbdIdxInternal(hist_key_indexs, hsv_key_indexs);
    //auto end = std::chrono::system_clock::now();
    //std::chrono::duration<double> elapsed_seconds = end - start;
    //std::cout<< "Tail time: " << elapsed_seconds.count() << "s\n";
    return sbd_indexs;
}

void DeepvacSB::reset() noexcept {
    hist_info_.clear();
    diff_vec_.clear();
    prev_hsv_pixels_.clear();
    hsv_diff_vec_.clear();
}

void DeepvacSB::appendFrameRGB(const unsigned char* rgb, const int width, const int height) noexcept {
    //auto start = std::chrono::system_clock::now();
    unsigned int result_rgb[DeepvacSB::channels_][DeepvacSB::hist_size_] = {0};
    width_ = width == 0 ? width_ : width;
    height_ = height ==0 ? height_ : height;
    //hist
    calcHistRGB(rgb, result_rgb);
    for(int i = 0; i < DeepvacSB::hist_size_; ++i) {
        for(int c=0; c<DeepvacSB::channels_; c++){
            hist_info_.push_back(result_rgb[c][i]);
        }
    }
    //hsv
    std::vector<unsigned char> hsv_pixels = rgb2Hsv();
    calcHSVDiff(hsv_pixels);

    //auto end = std::chrono::system_clock::now();
    //std::chrono::duration<double> elapsed_seconds = end - start;
    //std::cout<< "append time: " << elapsed_seconds.count() << "s\n";
}

void DeepvacSB::calcHistDiff() noexcept {
	const auto pixel_num = width_ * height_ * DeepvacSB::channels_;
	diff_vec_.clear();
	diff_vec_.push_back(0.0001);
	const auto frame_num = hist_info_.size() / DeepvacSB::hist_size_ / DeepvacSB::channels_;
	const auto temp = DeepvacSB::hist_size_*DeepvacSB::channels_;
	for(int i = 1; i < frame_num; ++i) {
		double diff_count = 0;
		for(int j = 0; j < DeepvacSB::hist_size_; ++j) {
			const auto prev_r = hist_info_[(i-1)*temp + j*3], next_r = hist_info_[i*temp + j*3];
			const auto prev_g = hist_info_[(i-1)*temp + j*3 + 1], next_g = hist_info_[i*temp + j*3 + 1];
			const auto prev_b = hist_info_[(i-1)*temp + j*3 + 2], next_b = hist_info_[i*temp + j*3 + 2];
			const auto v1 = next_r - prev_r;
			const auto v2 = next_g - prev_g;
			const auto v3 = next_b - prev_b;
			diff_count += (v1*v1) / (std::max(next_r, prev_r) + 1e-3);
			diff_count += (v2*v2) / (std::max(next_g, prev_g) + 1e-3);
			diff_count += (v3*v3) / (std::max(next_b, prev_b) + 1e-3);
		}

		diff_count /= pixel_num;
		diff_count *= 1;
		diff_count = std::max(0.0001, diff_count);
		diff_count = std::min(1., diff_count);
		diff_vec_.push_back(diff_count);
	}
}

std::vector<unsigned int> DeepvacSB::getSbdIdxInternal(const std::vector<unsigned int>& hist_index, const std::vector<unsigned int>& hsv_index) noexcept {
    std::vector<unsigned int> expand;
    for(int i = -3; i < 4; ++i) {
        for(int k = 0; k < hsv_index.size(); ++k) {
            int m = hsv_index[k]+i;
            expand.push_back( m >= 0 ? m:0);
        }
    }
    
    std::sort(expand.begin(), expand.end());
    std::vector<unsigned int> intersection;
    
    std::set_intersection(hist_index.begin(),hist_index.end(),
                          expand.begin(),expand.end(),
                          back_inserter(intersection));
    
    if(intersection.empty()) {
        return intersection;
    }

    std::vector<unsigned int> final{intersection[0]};
    for(int i = 1; i < intersection.size(); ++i) {
        if((intersection[i] - final.back()) > DeepvacSB::shot_diff_threshold_) {
            final.push_back(intersection[i]);
        }
    }

    return final;
}

std::vector<unsigned int> DeepvacSB::filterHistFrameIndex() noexcept {
	calcHistDiff();

    std::vector<unsigned int> indexs;
    const auto size = diff_vec_.size();
    float sum = 0;
    indexs.push_back(0);
    for (int i = 1; i < size; ++i) {
        sum += diff_vec_[i];
        auto avg = (sum) / (i);
        auto radio = fabs(diff_vec_[i] / (avg + (1e-3)));
        if (radio > DeepvacSB::hist_diff_threshold_) {
            indexs.push_back(i);
        }
    }
    return indexs;
}

std::vector<unsigned int> DeepvacSB::filterHsvFrameIndex() noexcept {
    std::vector<unsigned int> index_list{0};
    for(int i = 0; i < hsv_diff_vec_.size(); ++i) {
        if(hsv_diff_vec_[i] > DeepvacSB::hsv_diff_threshold_) {
            index_list.push_back(i);
        }
    }

    std::vector<unsigned int> key_index{0};
    for(int i = 1; i < index_list.size(); ++i) {
        if((index_list[i] - key_index[key_index.size()-1]) > DeepvacSB::shot_diff_threshold_) {
            key_index.push_back(index_list[i]);
        }
    }
    return key_index;
}

void DeepvacSB::downsampleImg(const unsigned char* colors) noexcept {
    const auto edge = std::max(width_, height_);
    auto factor_index = 0;
    for(; factor_index < DeepvacSB::level_; ++factor_index) {
    	if(edge >= resolution_ratio_level_[factor_index]) break;
    }

    const auto factor = DeepvacSB::downsample_factor_[factor_index];
    const auto new_width = width_ / factor;
    const auto new_height = height_ / factor;
    
    std::unique_ptr<unsigned char[]> tmp(new unsigned char[new_width * new_height * DeepvacSB::channels_]);
    rgb_colors_ = std::move(tmp);
    auto index = 0;
    for(int i = 0; i < height_; i += (factor)) {
    	for(int j = 0; j < width_; j += (factor)) {
            for(int k = 0; k < DeepvacSB::channels_; ++k) {
                rgb_colors_[index++] = colors[i*width_*DeepvacSB::channels_ + j*DeepvacSB::channels_ + k];
            }
    	}
    }
    //put in the end!!!
    width_ = new_width;
    height_ = new_height;
}

std::vector<unsigned char> DeepvacSB::rgb2Hsv() noexcept {
    std::vector<unsigned char> pixels;
    const auto size = height_*width_;
    auto index = 0;
    for(int i = 0; i < size; ++i) {
        int H, S, V;
        rgb2HsvSinglePixel(rgb_colors_[index+2], rgb_colors_[index+1], rgb_colors_[index], H, S, V);
        pixels.push_back(H);
        pixels.push_back(S);
        pixels.push_back(V);
        index += 3;
    }
    return pixels;
}

void DeepvacSB::rgb2HsvSinglePixel(const int R_I, const int G_I, const int B_I, int& H_I, int& S_I, int& V_I) noexcept {
    auto R = R_I / 255.f;
    auto G = G_I / 255.f;
    auto B = B_I / 255.f;
    float H, S, V;
    float min, max, delta,tmp;
    tmp = R>G?G:R;
    min = tmp>B?B:tmp;
    tmp = R>G?R:G;
    max = tmp>B?tmp:B;
    V_I = max * 255; // v
    V_I = std::min(V_I, 255);
    delta = max - min;
    if( max != 0 ) {
        S = delta / max; // s
        S_I = S * 255;
        S_I = std::min(S_I, 255);
    }else{
        S_I = 0;
        H_I = 0;
        return;
    }
    if (delta == 0){
        H_I = 0;
        return;
    }else if(R == max){
        if (G >= B){
            H = (G - B) / delta;
        }else{
            H = (G - B) / delta + 6.0;
        }
    }else if( G == max ){
        H = 2.0 + ( B - R ) / delta;
    }else if (B == max){
        H = 4.0 + ( R - G ) / delta;
    }
    H *= 60.0;
    H_I = H;
    H_I = std::min(H_I / 2, 179);
}

void DeepvacSB::calcHSVDiff(std::vector<unsigned char>& cur_hsv_pixels) noexcept {
    const auto size = prev_hsv_pixels_.size();
    if(!size) {
        hsv_diff_vec_.push_back(0);
        prev_hsv_pixels_.swap(cur_hsv_pixels);
        return;
    }

    int count = 0;
    int t = size / 4;
    int m = t * 4;
    for(int i = 0; i < m; i += 4) {
        count += std::abs((int)cur_hsv_pixels[i] - (int)prev_hsv_pixels_[i]);
        count += std::abs((int)cur_hsv_pixels[i+1] - (int)prev_hsv_pixels_[i+1]);
        count += std::abs((int)cur_hsv_pixels[i+2] - (int)prev_hsv_pixels_[i+2]);
        count += std::abs((int)cur_hsv_pixels[i+3] - (int)prev_hsv_pixels_[i+3]);
    }

    for(int i = m; m < size; ++m) {
        count += std::abs((int)cur_hsv_pixels[i] - (int)prev_hsv_pixels_[i]);
    }

    hsv_diff_vec_.push_back(count / float(size) * 1);
    prev_hsv_pixels_.swap(cur_hsv_pixels);
}
