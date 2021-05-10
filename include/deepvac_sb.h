/*
 * Copyright (c) 2020 DeepVAC
 * This file is part of libdeepvac, licensed under the GPLv3 (the "License")
 * You may not use this file except in compliance with the License.
 */

#pragma once
#include <algorithm>
#include <vector>
#include <iostream>
#include <memory>

class DeepvacSB {
public:
    DeepvacSB() = default;
    ~DeepvacSB() {
        reset();
    }

    void reset() noexcept;
    void appendFrameRGB(const unsigned char* rgb, const int width, const int height) noexcept;
    std::vector<unsigned int> getSbdIdx() noexcept;

private:
    static const int channels_{3};
    static const int hist_size_{256};
    static constexpr float hist_diff_threshold_{2.7f};

    static const int hsv_diff_threshold_{25};
    static const int shot_diff_threshold_{21};
    static const int level_{7};

    static const int resolution_ratio_level_[DeepvacSB::level_];
    static const int downsample_factor_[DeepvacSB::level_+1];

    int width_{0};
    int height_{0};
    int pixel_num_{0};
    int channel_pixel_num_{0};

    std::vector<float> hist_info_;
    std::vector<float> diff_vec_;

    std::vector<unsigned char> prev_hsv_pixels_;
    std::vector<float> hsv_diff_vec_;

    std::unique_ptr<unsigned char[]> rgb_colors_;


private:
    std::vector<unsigned int> getSbdIdxInternal(const std::vector<unsigned int>& hist_index, const std::vector<unsigned int>& hsv_index) noexcept;
    void downsampleImg(const unsigned char* colors) noexcept;

    //hist
    std::vector<unsigned int> filterHistFrameIndex() noexcept;
    void calcHistDiff() noexcept;

    //hsv
    std::vector<unsigned int> filterHsvFrameIndex() noexcept ;
    void calcHSVDiff(std::vector<unsigned char>& cur_hsv_pixels) noexcept;
    std::vector<unsigned char> rgb2Hsv() noexcept;
    void rgb2HsvSinglePixel(const int R_I, const int G_I, const int B_I, int& H_I, int& S_I, int& V_I) noexcept;

    void calcHistRGB(const unsigned char* rgb, unsigned int result_rgb[DeepvacSB::channels_][DeepvacSB::hist_size_]) {
        //must downsample 1st !!!
        downsampleImg(rgb);
        const auto pixel_num = width_ * height_;
        int m = pixel_num / 4;
        int size = m * 4;
        int i = 0;
        auto colors = rgb_colors_.get();
        for(; i < size; i += 4) {
            for(int g=0; g<4;g++){
                for(int c=0; c<DeepvacSB::channels_; c++){
                    ++result_rgb[c][*(colors++)];
                }
            }
        }

        for(; i < pixel_num; ++i) {
            for(int c=0; c<DeepvacSB::channels_; c++){
                ++result_rgb[c][*(colors++)];
            }
        }
    }
};
