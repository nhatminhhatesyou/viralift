# Auto-split from the original monolithic streamlit_app.py.
import streamlit as st


def _theme_overrides() -> str:
    """Follow the browser/OS color scheme instead of maintaining app theme state."""
    return """
    <style>
        :root {
            color-scheme: light dark;
            --vl-bg: #f4f7f3;
            --vl-bg-2: #e8efe8;
            --vl-surface: #fbfcf8;
            --vl-surface-strong: #ffffff;
            --vl-surface-muted: #eef4ee;
            --vl-border: #d6ded3;
            --vl-border-strong: #bac9b7;
            --vl-text: #17221d;
            --vl-muted: #637269;
            --vl-faint: #859188;
            --vl-accent: #1d6c63;
            --vl-accent-2: #174f49;
            --vl-accent-soft: #dceee9;
            --vl-warn: #9a6a1f;
            --vl-danger: #a63a3a;
            --vl-danger-soft: #f5dfdc;
            --vl-ok-soft: #dceee3;
            --vl-shadow: 0 22px 70px rgba(42, 70, 58, 0.12);
            --vl-app-bg:
                radial-gradient(circle at 6% 4%, rgba(29, 108, 99, 0.14), transparent 30rem),
                radial-gradient(circle at 86% 0%, rgba(71, 105, 86, 0.12), transparent 28rem),
                linear-gradient(180deg, var(--vl-bg), #fbfcf8 46%, #f6f8f4);
            --vl-sidebar-bg:
                linear-gradient(180deg, rgba(251, 252, 248, 0.96), rgba(236, 243, 235, 0.98)),
                var(--vl-surface);
            --vl-card-bg: linear-gradient(180deg, rgba(255, 255, 255, 0.96), rgba(247, 250, 246, 0.96));
            --vl-hero-bg:
                radial-gradient(circle at 95% 10%, rgba(29, 108, 99, 0.15), transparent 16rem),
                linear-gradient(135deg, rgba(255, 255, 255, 0.92), rgba(239, 246, 238, 0.88));
            --vl-soft-panel: rgba(255, 255, 255, 0.72);
            --vl-soft-panel-2: rgba(255, 255, 255, 0.74);
            --vl-upload-zone: rgba(238, 244, 238, 0.58);
            --vl-upload-zone-hover: rgba(220, 238, 233, 0.68);
            --vl-orb: rgba(29, 108, 99, 0.08);
        }

        @media (prefers-color-scheme: dark) {
            :root {
                --vl-bg: #07110f;
                --vl-bg-2: #0b1916;
                --vl-surface: #0d1b18;
                --vl-surface-strong: #10231f;
                --vl-surface-muted: #132d28;
                --vl-border: #24413b;
                --vl-border-strong: #38625a;
                --vl-text: #eef8f1;
                --vl-muted: #a1b4ab;
                --vl-faint: #6f877d;
                --vl-accent: #4fe0c6;
                --vl-accent-2: #21a894;
                --vl-accent-soft: rgba(79, 224, 198, 0.16);
                --vl-warn: #ffd37a;
                --vl-danger: #ff8f7f;
                --vl-danger-soft: rgba(255, 143, 127, 0.16);
                --vl-ok-soft: rgba(79, 224, 198, 0.14);
                --vl-shadow: 0 24px 80px rgba(0, 0, 0, 0.34);
                --vl-app-bg:
                    radial-gradient(circle at 8% 0%, rgba(79, 224, 198, 0.16), transparent 31rem),
                    radial-gradient(circle at 88% 3%, rgba(61, 126, 111, 0.18), transparent 30rem),
                    linear-gradient(180deg, #06100e, #091614 42%, #0b1110);
                --vl-sidebar-bg:
                    linear-gradient(180deg, rgba(12, 28, 24, 0.98), rgba(7, 17, 15, 0.98)),
                    var(--vl-surface);
                --vl-card-bg: linear-gradient(180deg, rgba(19, 42, 37, 0.94), rgba(10, 26, 23, 0.96));
                --vl-hero-bg:
                    radial-gradient(circle at 90% 7%, rgba(79, 224, 198, 0.18), transparent 18rem),
                    linear-gradient(135deg, rgba(18, 43, 38, 0.96), rgba(7, 18, 16, 0.93));
                --vl-soft-panel: rgba(16, 35, 31, 0.72);
                --vl-soft-panel-2: rgba(14, 31, 28, 0.82);
                --vl-upload-zone: rgba(17, 42, 37, 0.68);
                --vl-upload-zone-hover: rgba(27, 65, 58, 0.78);
                --vl-orb: rgba(79, 224, 198, 0.1);
            }
        }

        html, body, .stApp,
        [data-testid="stAppViewContainer"] {
            color-scheme: light dark;
            color: var(--vl-text) !important;
            background: var(--vl-app-bg) !important;
        }

        [data-testid="stHeader"] {
            background: transparent !important;
        }

        section[data-testid="stSidebar"],
        section[data-testid="stSidebar"] > div {
            color: var(--vl-text) !important;
            background: var(--vl-sidebar-bg) !important;
        }
    </style>
    """


def _inject_css():
    st.markdown(
        """
        <style>
            :root {
                --vl-bg: #f4f7f3;
                --vl-bg-2: #e8efe8;
                --vl-surface: #fbfcf8;
                --vl-surface-strong: #ffffff;
                --vl-surface-muted: #eef4ee;
                --vl-border: #d6ded3;
                --vl-border-strong: #bac9b7;
                --vl-text: #17221d;
                --vl-muted: #637269;
                --vl-faint: #859188;
                --vl-accent: #1d6c63;
                --vl-accent-2: #174f49;
                --vl-accent-soft: #dceee9;
                --vl-warn: #9a6a1f;
                --vl-danger: #a63a3a;
                --vl-danger-soft: #f5dfdc;
                --vl-ok-soft: #dceee3;
                --vl-shadow: 0 22px 70px rgba(42, 70, 58, 0.12);
            }

            html, body, [class*="css"], .stApp {
                font-family: "Geist", ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
                font-variant-numeric: tabular-nums;
            }

            .stApp {
                color: var(--vl-text);
                background: var(--vl-app-bg);
            }

            .main .block-container {
                max-width: 1320px;
                padding: 2.1rem 2.5rem 4rem;
            }

            section[data-testid="stSidebar"] {
                background: var(--vl-sidebar-bg);
                border-right: 1px solid var(--vl-border);
            }

            section[data-testid="stSidebar"] h1 {
                font-size: 1.28rem;
                letter-spacing: 0;
                font-weight: 800;
            }

            section[data-testid="stSidebar"] [data-testid="stMarkdownContainer"] p {
                color: var(--vl-muted);
            }

            h1, h2, h3 {
                color: var(--vl-text);
                letter-spacing: 0;
            }

            h1 {
                font-size: clamp(2.35rem, 4vw, 4.4rem);
                line-height: 0.96;
                margin-bottom: 0.65rem;
                font-weight: 800;
            }

            h2, h3 {
                font-weight: 750;
            }

            hr {
                margin: 1.45rem 0;
                border-color: var(--vl-border);
            }

            div[data-testid="stFileUploader"] {
                background: var(--vl-card-bg);
                border: 1px solid var(--vl-border);
                border-radius: 14px;
                padding: 1rem;
                box-shadow: 0 16px 45px rgba(45, 74, 63, 0.07);
            }

            div[data-testid="stFileUploader"] section {
                border: 1.4px dashed var(--vl-border-strong);
                border-radius: 12px;
                background: var(--vl-upload-zone);
                transition: border-color 160ms ease, background 160ms ease, transform 160ms ease;
            }

            div[data-testid="stFileUploader"] section:hover {
                border-color: var(--vl-accent);
                background: var(--vl-upload-zone-hover);
                transform: translateY(-1px);
            }

            div[data-testid="stFileUploader"] button,
            div[data-testid="stFileUploader"] button[kind],
            div[data-testid="stFileUploader"] [data-testid="stBaseButton-secondary"] {
                background: linear-gradient(180deg, var(--vl-surface-muted), var(--vl-surface-strong));
                border: 1px solid var(--vl-border-strong);
                color: var(--vl-text);
                border-radius: 10px;
                box-shadow: none;
                font-weight: 750;
            }

            div[data-testid="stFileUploader"] button:hover,
            div[data-testid="stFileUploader"] [data-testid="stBaseButton-secondary"]:hover {
                background: var(--vl-accent-soft);
                border-color: var(--vl-accent);
                color: var(--vl-text);
                box-shadow: 0 10px 22px rgba(79, 224, 198, 0.13);
            }

            div[data-testid="stFileUploader"] button *,
            div[data-testid="stFileUploader"] [data-testid="stBaseButton-secondary"] * {
                color: var(--vl-text);
                fill: var(--vl-text);
                stroke: var(--vl-text);
            }

            div[data-testid="stFileUploader"] small,
            div[data-testid="stFileUploader"] span,
            div[data-testid="stFileUploader"] p {
                color: var(--vl-muted);
            }

            label, .stMarkdown p, .stCaptionContainer {
                color: var(--vl-muted);
            }

            div[data-testid="stMetric"] {
                background: var(--vl-card-bg), var(--vl-surface);
                border: 1px solid var(--vl-border);
                border-radius: 14px;
                padding: 1rem 1.05rem;
                box-shadow: var(--vl-shadow);
                min-height: 7rem;
            }

            div[data-testid="stMetricLabel"] p {
                color: var(--vl-muted);
                font-size: 0.76rem;
                font-weight: 700;
            }

            div[data-testid="stMetricValue"] {
                color: var(--vl-text);
                font-weight: 800;
            }

            div[data-testid="stMetricDelta"] {
                color: var(--vl-danger);
            }

            .vl-kicker {
                color: var(--vl-accent);
                font-size: 0.72rem;
                font-weight: 800;
                margin-bottom: 0.55rem;
                text-transform: uppercase;
                letter-spacing: 0.12em;
            }

            .vl-help {
                color: var(--vl-muted);
                font-size: 1.02rem;
                max-width: 820px;
                line-height: 1.62;
                margin-bottom: 1.1rem;
            }

            .vl-hero {
                position: relative;
                overflow: hidden;
                border: 1px solid var(--vl-border);
                border-radius: 18px;
                padding: 1.25rem 1.35rem 1.45rem;
                margin-bottom: 1.35rem;
                background: var(--vl-hero-bg);
                box-shadow: var(--vl-shadow);
            }

            .vl-hero::after {
                content: "";
                position: absolute;
                inset: auto -6rem -7rem auto;
                width: 18rem;
                height: 18rem;
                border-radius: 999px;
                background: var(--vl-orb);
                pointer-events: none;
            }

            .vl-stage-rail {
                position: relative;
                margin: 1.05rem 0 0.1rem;
                max-width: 54rem;
            }

            .vl-stage-meta {
                color: var(--vl-accent);
                font-size: 0.72rem;
                font-weight: 800;
                letter-spacing: 0.09em;
                text-transform: uppercase;
                margin-bottom: 0.7rem;
            }

            .vl-stage-track {
                position: relative;
                display: grid;
                grid-template-columns: repeat(5, minmax(0, 1fr));
                gap: 0;
                list-style: none;
                padding: 0;
                margin: 0;
            }

            /* Connector: a muted base line plus an accent fill up to the
               active step (width driven by --vl-progress set inline). Both run
               between the first and last dot centres (10% .. 90%). */
            .vl-stage-track::before,
            .vl-stage-track::after {
                content: "";
                position: absolute;
                left: 10%;
                top: 0.88rem;
                height: 2px;
                border-radius: 999px;
            }

            .vl-stage-track::before {
                right: 10%;
                background: var(--vl-border-strong);
                opacity: 0.55;
            }

            .vl-stage-track::after {
                width: calc(var(--vl-progress, 10%) - 10%);
                background: linear-gradient(
                    90deg,
                    color-mix(in srgb, var(--vl-accent) 60%, transparent),
                    var(--vl-accent)
                );
                box-shadow: 0 0 10px color-mix(in srgb, var(--vl-accent) 45%, transparent);
                transition: width 360ms cubic-bezier(0.22, 0.61, 0.36, 1);
            }

            .vl-step {
                position: relative;
                z-index: 1;
                min-width: 0;
                text-align: center;
                color: var(--vl-faint);
            }

            .vl-step-dot {
                position: relative;
                display: inline-grid;
                place-items: center;
                width: 1.8rem;
                height: 1.8rem;
                border-radius: 999px;
                border: 1.5px solid var(--vl-border-strong);
                background: var(--vl-surface);
                color: var(--vl-muted);
                font-size: 0.76rem;
                font-weight: 850;
                transition: transform 200ms ease, background 200ms ease,
                            border-color 200ms ease;
            }

            .vl-step-label {
                display: block;
                margin: 0.5rem auto 0;
                max-width: 8rem;
                color: inherit;
                font-size: 0.74rem;
                font-weight: 700;
                line-height: 1.2;
                white-space: nowrap;
                overflow: hidden;
                text-overflow: ellipsis;
            }

            /* completed: filled accent dot with a check mark */
            .vl-step-complete {
                color: var(--vl-muted);
            }

            .vl-step-complete .vl-step-dot {
                border-color: var(--vl-accent);
                background: var(--vl-accent);
                color: transparent;
            }

            .vl-step-complete .vl-step-dot::after {
                content: "✓";
                position: absolute;
                inset: 0;
                display: grid;
                place-items: center;
                color: var(--vl-bg);
                font-size: 0.9rem;
                font-weight: 900;
                line-height: 1;
            }

            /* active: bold accent dot with a glow ring */
            .vl-step-active {
                color: var(--vl-text);
            }

            .vl-step-active .vl-step-dot {
                border-color: var(--vl-accent);
                background: var(--vl-accent);
                color: var(--vl-bg);
                transform: scale(1.14);
                box-shadow:
                    0 0 0 5px color-mix(in srgb, var(--vl-accent) 20%, transparent),
                    0 8px 22px color-mix(in srgb, var(--vl-accent) 30%, transparent);
            }

            .vl-step-active .vl-step-label {
                color: var(--vl-text);
                font-weight: 850;
            }

            .vl-panel {
                background: var(--vl-card-bg);
                border: 1px solid var(--vl-border);
                border-radius: 14px;
                padding: 1.05rem 1.15rem;
                margin: 0.6rem 0 1.1rem;
                box-shadow: 0 14px 40px rgba(45, 74, 63, 0.06);
            }

            .vl-panel strong {
                color: var(--vl-text);
            }

            .vl-pill {
                display: inline-block;
                border: 1px solid var(--vl-border);
                border-radius: 999px;
                padding: 0.22rem 0.55rem;
                margin: 0.1rem;
                background: var(--vl-surface-muted);
                color: var(--vl-text);
                font-size: 0.82rem;
                font-family: "Geist Mono", ui-monospace, SFMono-Regular, Menlo, monospace;
            }

            .vl-context-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(8.5rem, 1fr));
                gap: 0.6rem;
                margin: 1rem 0 1.25rem;
            }

            .vl-context {
                background: var(--vl-soft-panel);
                border: 1px solid var(--vl-border);
                border-radius: 12px;
                padding: 0.85rem 0.95rem;
            }

            .vl-context-label {
                color: var(--vl-faint);
                font-size: 0.7rem;
                font-weight: 800;
                letter-spacing: 0.1em;
                text-transform: uppercase;
            }

            .vl-context-value {
                color: var(--vl-text);
                font-family: "Geist Mono", ui-monospace, SFMono-Regular, Menlo, monospace;
                font-size: 0.88rem;
                font-weight: 700;
                margin-top: 0.28rem;
                overflow-wrap: anywhere;
            }

            .vl-sidebar-card {
                border: 1px solid var(--vl-border);
                border-radius: 14px;
                padding: 0.9rem;
                margin: 0.75rem 0;
                background: var(--vl-soft-panel);
            }

            .vl-sidebar-label {
                color: var(--vl-faint);
                font-size: 0.68rem;
                font-weight: 800;
                letter-spacing: 0.1em;
                text-transform: uppercase;
            }

            .vl-sidebar-value {
                color: var(--vl-text);
                font-family: "Geist Mono", ui-monospace, SFMono-Regular, Menlo, monospace;
                font-size: 0.82rem;
                font-weight: 700;
                margin-top: 0.22rem;
                overflow-wrap: anywhere;
            }

            .stButton > button,
            .stDownloadButton > button {
                border-radius: 10px;
                border: 1px solid var(--vl-border-strong);
                color: var(--vl-text);
                font-weight: 750;
                transition: transform 160ms ease, border-color 160ms ease, background 160ms ease, box-shadow 160ms ease;
            }

            .stButton > button *,
            .stDownloadButton > button * {
                color: inherit !important;
                -webkit-text-fill-color: currentColor !important;
            }

            .stButton > button:hover,
            .stDownloadButton > button:hover {
                transform: translateY(-1px);
                border-color: var(--vl-accent);
                box-shadow: 0 12px 24px rgba(29, 108, 99, 0.12);
            }

            .stButton > button:active,
            .stDownloadButton > button:active {
                transform: translateY(0);
            }

            button[kind="primary"],
            .stButton > button[kind="primary"] {
                border-color: var(--vl-accent-2);
                background: linear-gradient(180deg, var(--vl-accent), var(--vl-accent-2));
                color: #ffffff !important;
                -webkit-text-fill-color: #ffffff !important;
            }

            button[kind="primary"] *,
            .stButton > button[kind="primary"] * {
                color: #ffffff !important;
                -webkit-text-fill-color: #ffffff !important;
            }

            .st-key-delete_virus_entry button {
                border-color: var(--vl-danger) !important;
                background: linear-gradient(180deg, var(--vl-danger), #7f2525) !important;
                color: #ffffff !important;
            }

            .st-key-delete_virus_entry button:hover {
                border-color: var(--vl-danger) !important;
                box-shadow: 0 14px 28px rgba(166, 58, 58, 0.2) !important;
            }

            div[data-testid="stExpander"] {
                border: 1px solid var(--vl-border);
                border-radius: 13px;
                background: var(--vl-soft-panel-2);
                box-shadow: 0 12px 34px rgba(45, 74, 63, 0.05);
            }

            div[data-testid="stExpander"] summary {
                font-weight: 750;
                color: var(--vl-text);
            }

            div[data-testid="stDataFrame"] {
                border: 1px solid var(--vl-border);
                border-radius: 13px;
                overflow: hidden;
                box-shadow: 0 16px 40px rgba(45, 74, 63, 0.06);
            }

            /* Spread the tabs evenly across the full width like a segmented bar
               instead of clustering on the left. */
            div[data-baseweb="tab-list"] {
                display: flex;
                width: 100%;
                gap: 0.4rem;
                border-bottom: 1px solid var(--vl-border);
                margin-bottom: 1.4rem;
            }

            button[data-baseweb="tab"] {
                flex: 1 1 0;
                justify-content: center;
                text-align: center;
                border-radius: 10px 10px 0 0;
                font-weight: 750;
                padding: 0.4rem 0.5rem 0.7rem;
            }

            /* breathing room between the tab bar and the panel content */
            div[data-baseweb="tab-panel"] {
                padding-top: 0.4rem;
            }

            div[data-testid="stAlert"] {
                border-radius: 12px;
                border: 1px solid var(--vl-border);
            }

            .stProgress > div > div > div {
                background-color: var(--vl-accent);
            }

            input, textarea, div[data-baseweb="select"] > div {
                border-radius: 10px;
            }

            div[role="radiogroup"] {
                gap: 0.35rem;
            }

            div[role="radiogroup"] label {
                border: 1px solid var(--vl-border);
                border-radius: 999px;
                padding: 0.25rem 0.65rem;
                background: var(--vl-soft-panel);
                margin-right: 0.25rem;
            }

            div[role="radiogroup"] label:has(input:checked) {
                border-color: var(--vl-accent);
                background: var(--vl-accent-soft);
                color: var(--vl-text);
            }

            div[data-baseweb="textarea"] > div,
            div[data-baseweb="input"] > div,
            div[data-baseweb="select"] > div {
                background: var(--vl-soft-panel);
                border-color: var(--vl-border);
                color: var(--vl-text);
            }

            textarea,
            textarea:focus,
            textarea:disabled,
            div[data-baseweb="textarea"] textarea,
            div[data-baseweb="textarea"] textarea:focus {
                background: transparent !important;
                color: var(--vl-text) !important;
                -webkit-text-fill-color: var(--vl-text) !important;
                caret-color: var(--vl-accent) !important;
            }

            div[data-baseweb="textarea"] {
                background: var(--vl-soft-panel) !important;
                border-radius: 10px;
            }

            div[data-baseweb="input"],
            div[data-baseweb="textarea"],
            div[data-baseweb="select"],
            div[data-baseweb="input"] *,
            div[data-baseweb="textarea"] *,
            div[data-baseweb="select"] * {
                color: var(--vl-text);
                caret-color: var(--vl-accent);
            }

            div[data-baseweb="input"] input,
            div[data-baseweb="textarea"] textarea {
                color: var(--vl-text);
                -webkit-text-fill-color: var(--vl-text);
            }

            div[data-baseweb="input"] input::placeholder,
            div[data-baseweb="textarea"] textarea::placeholder {
                color: var(--vl-faint);
                -webkit-text-fill-color: var(--vl-faint);
            }

            div[data-baseweb="tab-list"] button,
            div[data-baseweb="tab-list"] button * {
                color: var(--vl-muted);
            }

            div[data-baseweb="tab-list"] button[aria-selected="true"],
            div[data-baseweb="tab-list"] button[aria-selected="true"] * {
                color: var(--vl-danger);
            }

            div[data-baseweb="popover"] {
                background: var(--vl-surface-strong);
                color: var(--vl-text);
            }

            div[data-testid="stAlert"] *,
            div[data-testid="stNotification"] * {
                color: var(--vl-text);
            }

            @media (max-width: 760px) {
                .main .block-container {
                    padding: 1.15rem 1rem 2rem;
                }

                .vl-context-grid {
                    grid-template-columns: 1fr;
                }

                .vl-stage-track {
                    grid-template-columns: repeat(3, minmax(0, 1fr));
                    row-gap: 0.8rem;
                }

                .vl-stage-track::before,
                .vl-stage-track::after {
                    display: none;
                }

                h1 {
                    font-size: 2.35rem;
                }
            }
        </style>
        """,
        unsafe_allow_html=True,
    )
    st.markdown(_theme_overrides(), unsafe_allow_html=True)
