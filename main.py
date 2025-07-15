import streamlit as st
import numpy as np
from astropy.io import fits
from PIL import Image
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from datetime import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# --- Streamlit 앱 페이지 설정 ---
st.set_page_config(page_title="천문 이미지 분석기", layout="wide")
st.title("🔭 천문 이미지 처리 앱")

# --- 파일 업로더 ---
uploaded_file = st.file_uploader(
    "분석할 FITS 파일을 선택하세요.",
    type=['fits', 'fit', 'fz']
)

# --- 서울 위치 설정 (고정값) ---
seoul_location = EarthLocation(lat=37.5665, lon=126.9780, height=50)  # 서울 위도/경도/고도

# --- 현재 시간 (UTC 기준) ---
now = datetime.utcnow()
now_astropy = Time(now)

# --- 파일이 업로드되면 실행될 로직 ---
if uploaded_file:
    try:
        with fits.open(uploaded_file) as hdul:
            image_hdu = None
            for hdu in hdul:
                if hdu.data is not None and hdu.is_image:
                    image_hdu = hdu
                    break

            if image_hdu is None:
                st.error("파일에서 유효한 이미지 데이터를 찾을 수 없습니다.")
            else:
                header = image_hdu.header
                data = image_hdu.data
                data = np.nan_to_num(data)

                # --- FITS 헤더에서 거리 정보 추출 ---
                distance_keys = ['DISTANCE', 'DIST', 'PARSEC']
                header_distance = None
                for key in distance_keys:
                    if key in header:
                        try:
                            header_distance = float(header[key])
                            if header_distance <= 0:
                                header_distance = None
                                st.sidebar.warning(f"FITS 헤더의 {key} 값({header[key]})이 유효하지 않습니다. 기본값(100 파섹)을 사용합니다.")
                            break
                        except (ValueError, TypeError):
                            st.sidebar.warning(f"FITS 헤더의 {key} 값({header[key]})을 숫자로 변환할 수 없습니다. 기본값(100 파섹)을 사용합니다.")
                            header_distance = None

                # 기본 거리 설정
                default_distance = header_distance if header_distance is not None else 100.0

                st.success(f"**'{uploaded_file.name}'** 파일을 성공적으로 처리했습니다.")
                col1, col2 = st.columns(2)

                with col1:
                    st.header("이미지 정보")
                    st.text(f"크기: {data.shape[1]} x {data.shape[0]} 픽셀")
                    if 'OBJECT' in header:
                        st.text(f"관측 대상: {header['OBJECT']}")
                    if 'EXPTIME' in header:
                        st.text(f"노출 시간: {header['EXPTIME']} 초")

                    st.header("물리량")
                    mean_brightness = np.mean(data)
                    st.metric(label="이미지 전체 평균 밝기", value=f"{mean_brightness:.2f}")

                with col2:
                    st.header("이미지 미리보기")
                    if data.max() == data.min():
                        norm_data = np.zeros(data.shape, dtype=np.uint8)
                    else:
                        scale_min = np.percentile(data, 5)
                        scale_max = np.percentile(data, 99.5)
                        data_clipped = np.clip(data, scale_min, scale_max)
                        norm_data = (255 * (data_clipped - scale_min) / (scale_max - scale_min)).astype(np.uint8)

                    img = Image.fromarray(norm_data)
                    st.image(img, caption="업로드된 FITS 이미지", use_container_width=True)

                # --- 3D 공간 좌표 시각화 ---
               st.header("3D Spatial Coordinates")
if 'RA' in header and 'DEC' in header:
    try:
        # Distance input
        distance = st.sidebar.number_input(
            "Distance to Object (parsecs)",
            min_value=0.0,
            value=default_distance,
            step=10.0,
            help="Distance read from FITS header or default value of 100 parsecs."
        )
        if header_distance is not None:
            st.sidebar.info(f"Distance of {header_distance:.2f} parsecs read from FITS header.")
        if distance <= 0:
            st.sidebar.warning("Distance must be positive.")
            st.write("Please enter a valid distance to display the 3D plot.")
        else:
            target_coord = SkyCoord(
                ra=header['RA'], dec=header['DEC'], distance=distance,
                unit=('hourangle', 'deg', 'parsec')
            )
            # 3D coordinates calculation
            x, y, z = target_coord.cartesian.xyz.value  # Parsecs

            # 3D plot creation
            fig = plt.figure(figsize=(8, 8))
            ax = fig.add_subplot(111, projection='3d')
            # Sun (origin)
            ax.scatter([0], [0], [0], color='yellow', s=100, label='Sun')
            # Celestial object
            ax.scatter([x], [y], [z], color='blue', s=50, label='Object')
            ax.set_xlabel('X (parsecs)')
            ax.set_ylabel('Y (parsecs)')
            ax.set_zlabel('Z (parsecs)')
            ax.set_title('3D Position of Celestial Object')
            ax.legend()
            ax.grid(True)
            # Optional: Dynamic axis limits for better visualization
            max_range = max(abs(x), abs(y), abs(z), distance) * 1.2
            ax.set_xlim([-max_range, max_range])
            ax.set_ylim([-max_range, max_range])
            ax.set_zlim([-max_range, max_range])
            st.pyplot(fig)
    except Exception as e:
        st.write(f"Failed to create 3D plot: {e}")
else:
    st.write("RA/DEC information not found in FITS header. Unable to display 3D plot.")

                # --- 사이드바: 현재 천체 위치 계산 및 적경/적위/거리 표시 ---
                st.sidebar.header("🧭 현재 천체 위치 (서울 기준)")
                if 'RA' in header and 'DEC' in header:
                    try:
                        distance = st.sidebar.number_input(
                            "천체까지의 거리 (파섹)",
                            min_value=0.0,
                            value=default_distance,
                            step=10.0,
                            key="distance_input",
                            help="FITS 헤더에서 거리 정보를 읽었거나 기본값 100 파섹을 사용합니다."
                        )
                        if distance <= 0:
                            target_coord = SkyCoord(
                                ra=header['RA'], dec=header['DEC'], unit=('hourangle', 'deg')
                            )
                        else:
                            target_coord = SkyCoord(
                                ra=header['RA'], dec=header['DEC'], distance=distance,
                                unit=('hourangle', 'deg', 'parsec')
                            )
                        altaz = target_coord.transform_to(AltAz(obstime=now_astropy, location=seoul_location))
                        altitude = altaz.alt.degree
                        azimuth = altaz.az.degree

                        # 적경과 적위 표시
                        ra_hms = target_coord.ra.to_string(unit='hour', sep=':', precision=2)
                        dec_dms = target_coord.dec.to_string(unit='deg', sep=':', precision=2)
                        st.sidebar.metric("적경 (RA)", ra_hms)
                        st.sidebar.metric("적위 (DEC)", dec_dms)
                        # 거리 표시
                        st.sidebar.metric("거리 (parsec)", f"{distance:.2f}" if distance > 0 else "N/A")
                        st.sidebar.metric("고도 (°)", f"{altitude:.2f}")
                        st.sidebar.metric("방위각 (°)", f"{azimuth:.2f}")
                    except Exception as e:
                        st.sidebar.warning(f"천체 위치 계산 실패: {e}")
                else:
                    st.sidebar.info("FITS 헤더에 RA/DEC 정보가 없습니다.")

    except Exception as e:
        st.error(f"파일 처리 중 오류가 발생했습니다: {e}")
        st.warning("파일이 손상되었거나 유효한 FITS 형식이 아닐 수 있습니다.")
else:
    st.info("시작하려면 FITS 파일을 업로드해주세요.")

# --- 💬 댓글 기능 (세션 기반) ---
st.divider()
st.header("💬 의견 남기기")

if "comments" not in st.session_state:
    st.session_state.comments = []

with st.form(key="comment_form"):
    name = st.text_input("이름을 입력하세요", key="name_input")
    comment = st.text_area("댓글을 입력하세요", key="comment_input")
    submitted = st.form_submit_button("댓글 남기기")

    if submitted:
        if name.strip() and comment.strip():
            st.session_state.comments.append((name.strip(), comment.strip()))
            st.success("댓글이 저장되었습니다.")
        else:
            st.warning("이름과 댓글을 모두 입력해주세요.")

if st.session_state.comments:
    st.subheader("📋 전체 댓글")
    for i, (n, c) in enumerate(reversed(st.session_state.comments), 1):
        st.markdown(f"**{i}. {n}**: {c}")
else:
    st.info("아직 댓글이 없습니다. 첫 댓글을 남겨보세요!")
